//
// Created by erolsson on 2019-01-20.
//

#include "stone_material_contact.h"

#include <random>

#include "../materials/stone_material.h"

DEM::StoneMaterialContact::StoneMaterialContact(ParticleType* particle1, ParticleType* particle2,
                                                std::chrono::duration<double>) :
                                                p1_{particle1}, p2_{particle2}
{
    auto mat1 = dynamic_cast<const StoneMaterial*>(particle1->get_material());
    auto mat2 = dynamic_cast<const StoneMaterial*>(particle2->get_material());

    double E1 = mat1->E;
    double E2 = mat2->E;
    double v1 = mat1->nu;
    double v2 = mat2->nu;
    double E0 = 1/((1 - v1*v1)/E1 + (1 - v2*v2)/E2);

    double R1 = particle1->get_radius();
    double R2 = particle2->get_radius();

    R0_ = 1/(1/R1 + 1/R2);

    kp_ = E0*4./3*sqrt(R0_);           //Loading stiffness
    ke_ = E0*4./3*sqrt(R0_)/0.95;
    kl_ = E0*4./3*sqrt(R0_)/0.95;

    hs_ = (mat1->hs + mat2->hs)/2;
    Fs_ = (mat1->Fs + mat2->Fs)/2/sqrt(0.00625/R0_);

    hlinear_ = (mat1->plastic_linear_depth + mat2->plastic_linear_depth)/2*R0_/6.25e-3 + hs_;

    h1_ = hs_ - pow(Fs_/kp_, 2./3);

    unloading_exp_ = (mat1->unloading_exponent + mat2->unloading_exponent)/2;
    ku_ = E0*4./3*pow(R0_, 2-unloading_exp_)/0.95;
    mu_ = (mat1->mu + mat2->mu)/2;
    old_mu_ = (mat1->mu + mat2->mu)/2;

    double G1 = E1/2/(1+v1);
    double G2 = E2/2/(1+v2);
    kT_ = 8/((2-v1)/G1 + (2-v2)/G2);

    id2_ = p2_->get_id();
    assign_fracture_strengths();
}

DEM::StoneMaterialContact::StoneMaterialContact(ParticleType* particle1, SurfaceType* surface,
                                                std::chrono::duration<double>) : p1_{particle1}, p2_{nullptr}
{
    auto mat1 = dynamic_cast<const StoneMaterial*>(particle1->get_material());

    double E1 = mat1->E;
    double v1 = mat1->nu;
    double E0 = 1/((1 - v1*v1)/E1);

    double R1 = particle1->get_radius();
    R0_ = R1;

    kp_ = E0*4./3*sqrt(R0_);           //Loading stiffness
    ke_ = E0*4./3*sqrt(R0_)/0.95;
    kl_ = E0*4./3*sqrt(R0_)/0.95;

    hs_ = mat1->hs;
    Fs_ = mat1->Fs/sqrt(0.00625/R0_);
    h1_ = hs_ - pow(Fs_/kp_, 2./3);

    hlinear_ = mat1->plastic_linear_depth*R0_/6.25e-3 + hs_;

    unloading_exp_ = mat1->unloading_exponent;
    ku_ = E0*4./3*pow(R0_, 2-unloading_exp_)/0.95;           //Unloading stiffness

    mu_ = mat1->mu_wall;
    old_mu_ = mat1->mu_wall;

    double G1 = E1/2/(1+v1);
    kT_ = 8/((2-v1)/G1);
    assign_fracture_strengths();

    id2_ = surface->get_id();
}

DEM::StoneMaterialContact::StoneMaterialContact(ParticleType* p1, ParticleType* p2, SurfaceType*,
    std::chrono::duration<double>,  const ParameterMap &parameters) :
    p1_(p1),
    p2_(p2),
    FN_(parameters.get_parameter("FN")),
    a_(parameters.get_parameter("a")),
    R0_(parameters.get_parameter("R0")),
    Fs_(parameters.get_parameter("Fs")),
    hs_(parameters.get_parameter("hs")),
    hlinear_(parameters.get_parameter("hlinear")),
    h_(parameters.get_parameter("h")),
    hmax_(parameters.get_parameter("hmax")),
    Fmax_(parameters.get_parameter("Fmax")),
    hp_(parameters.get_parameter("hp")),
    hl_(parameters.get_parameter("hl")),
    h1_(parameters.get_parameter("h1")),
    FT_(parameters.get_vec3("FT")),
    old_dT_(parameters.get_vec3("old_dT")),
    load_var_tangential_(parameters.get_vec3("load_var_tangential")),
    load_var_normal_(parameters.get_parameter("load_var_normal")),
    turning_point_(parameters.get_parameter<int>("turning_point")),
    FT0({parameters.get_vec3("FTO_1"), parameters.get_vec3("FTO_2")}),
    kp_(parameters.get_parameter("kp")),
    ke_(parameters.get_parameter("ke")),
    kl_(parameters.get_parameter("kl")),
    ku_(parameters.get_parameter("ku")),
    unloading_exp_(parameters.get_parameter("unloading_exp")),
    mu_(parameters.get_parameter("mu")),
    old_mu_(parameters.get_parameter("old_mu")),
    kT_(parameters.get_parameter("kT")),
    mu_r_(parameters.get_parameter("mu_r")),
    kr_(parameters.get_parameter("kr")),
    Tr_(parameters.get_vec3("Tr")),
    fracture_strength_p1_(parameters.get_parameter("fracture_strength_p1")),
    fracture_strength_p2_(parameters.get_parameter("fracture_strength_p2")),
    id2_(parameters.get_parameter<std::size_t>("id2"))
{

}

void DEM::StoneMaterialContact::update(double h, const DEM::Vec3& dt, const Vec3& rot, const DEM::Vec3& normal)
{
    double new_FN = update_normal_force(h - h_);
    double d_mu_FN = new_FN*mu_ - FN_*old_mu_;
    FN_ = new_FN;
    update_tangential_force(dt, normal, d_mu_FN);
    update_tangential_resistance(rot);
    check_fracture(normal);
}

double DEM::StoneMaterialContact::update_normal_force(double dh) {
    h_ += dh;
    double F = 0;
    if (h_ > 0) {
        a_ = sqrt(h_*R0_);
        if (h_ >= hmax_) {
            hmax_ = h_;
            if (h_ < hs_) {
                F = Fs_/hs_*h_;
            }
            else if (h_ < hlinear_){
                F = kp_*pow(h_ - h1_, 1.5);
            }
            else {
                double Flin = kp_*pow(hlinear_ - h1_, 1.5);
                double klin = 1.5*kp_*sqrt(hlinear_ - h1_);
                F = Flin + klin*(h_ - hlinear_);
            }
            Fmax_ = F;
            hp_ = (h_ - pow(F/ke_, 1./1.5)); // Plastic indentation depth
            if (hp_ < 0) {
                hp_ = hmax_/2;
            }
            ku_ = F/pow(h_-hp_, unloading_exp_);
        }

        else if (dh > 0) {
            if (hmax_ < hs_) {
                F = kl_*(h_ - hl_);
                hp_ = (h_ - hl_)/2;
            }
            else if (h_ > hl_) {
                F = kl_*pow(h_ - hl_, 1.5);
            }
            if (h_ > hp_) {
                ku_ = F/pow(h_ - hp_, unloading_exp_);
            }
        }
        else {
            if (h_ > hp_) {
                F = ku_*pow(h_ - hp_, unloading_exp_);
            }
            if (hmax_ < hs_) {
                kl_ = (Fmax_ - F)/(hmax_ - h_);
                hl_ = hmax_ - Fmax_/kl_;
            }
            else {
                double A = pow(F/Fmax_, 1./1.5);
                if (A > 1.) {
                    A = 1 - 1e-9;
                }
                hl_ = (h_ - A*hmax_)/(1 - A);
                kl_ = Fmax_/pow(hmax_ - hl_, 1.5);
            }
        }
    }
    else {
        a_ = 0;
        kl_ = ke_;
        hl_ = hp_;
        hmax_ = 0;
        hp_ = 0;
    }
    return F;
}

void DEM::StoneMaterialContact::update_tangential_force(const Vec3& dt, const Vec3& normal, double d_mu_FN) {
    /*
     * Friction model according to simplified Mindlin in
     * "An investigation of the comparative behaviour of alternative contact force models during elastic collisions"
     * All equation numbers refer to that paper
    */

    if (mu_*FN_ > 0) {  // Only care about friction if we have a normal force
        if (dt.length() > 1e-40) {
            // Project previous contact force on the contact plane
            FT_ -= dot_product(FT_, normal)*normal;
            // We have a change in loading direction, from "loading" to unloading
            if (dot_product(dt, old_dT_) < 0.) {
                turning_point_ *= -1;
                FT0[(turning_point_ + 1)/2] = FT_;
            }

            // Update the turning points with changes in normal force
            for (auto &F0: FT0) {
                if (F0.length() > 1e-40) {
                    // Equation (23) for the 2D case F** is negative and a - sign is used
                    // This is accounted by the normal instead of F0
                    F0 += d_mu_FN*F0.normal();
                }
            }
            int sgn = 1;
            double q = 0;
            if (d_mu_FN <= 0 || load_var_tangential_.length() > load_var_normal_) {
                if (dot_product(dt, old_dT_) >= 0. && FT0[0].is_zero()) {     // Loading of the contact
                    q = 1 - (FT_.length() + d_mu_FN)/(mu_*FN_);
                }
                else if (turning_point_ < 0) {                    // Unloading
                    q = 1 - ((- FT0[0] + FT_).length() + 2*d_mu_FN)/(2*mu_*FN_);
                    sgn = -1;
                }
                else {                                                        // Reloading
                    q = 1 - ((- FT_ + FT0[1]).length() + 2*d_mu_FN)/(2*mu_*FN_);
                }
                if (d_mu_FN < 0) {
                    load_var_tangential_.set_zero();
                    load_var_normal_ = 0;
                }
            }
            else {
                q = 1;
                load_var_tangential_ += kT_*a_*dt;
                load_var_normal_ += d_mu_FN;
            }
            if (q > 0) {
                q = pow(q, 1./3);
            }
            else {
                q = 0;
            }

            double kt = kT_*a_*q + sgn*d_mu_FN*(1 - q)/dt.length();
            FT_ -= kt*dt;
            old_dT_ = dt;
            if (FT_.length() > mu_*FN_)
                FT_ = mu_*FN_*FT_.normalize();
        }
    }
    else  {
        // Reset all state variables
        FT_.set_zero();
        for (auto& F0: FT0) {
            F0.set_zero();
        }
        turning_point_ = 1;
        old_dT_.set_zero();
        load_var_tangential_.set_zero();
        load_var_normal_ = 0;
    }
}

void DEM::StoneMaterialContact::update_tangential_resistance(const DEM::Vec3&) {
    // rot.length();
}

void DEM::StoneMaterialContact::assign_fracture_strengths() {
    std::default_random_engine random_engine(std::chrono::system_clock::now().time_since_epoch().count());
    auto mat1 = dynamic_cast<const StoneMaterial*>(p1_->get_material());
    std::weibull_distribution distr1(mat1->weibull_exponent,
                                             mat1->weibull_fracture_stress/pow(R0_, 3)*mat1->weibull_ref_volume);
    fracture_strength_p1_ = distr1(random_engine)*R0_*R0_;
    if(p2_ != nullptr) {
        auto mat2 = dynamic_cast<const StoneMaterial *>(p2_->get_material());
        std::weibull_distribution distr2(mat2->weibull_exponent,
                                         mat2->weibull_fracture_stress/pow(R0_, 3)*mat2->weibull_ref_volume);
        fracture_strength_p2_ = distr2(random_engine)*R0_*R0_;
    }
}

void DEM::StoneMaterialContact::check_fracture(const Vec3& normal) {
    if (FN_ > fracture_strength_p1_) {
        Vec3 contact_position = p1_->get_position() - p1_->get_radius()*normal;
        p1_->fracture_particle(contact_position, FN_, id2_, normal);
    }
    if (p2_ != nullptr && FN_ > fracture_strength_p2_) {
        Vec3 contact_position = p2_->get_position() + p2_->get_radius()*normal;
        p2_->fracture_particle(contact_position, FN_, p1_->get_id(), normal);
    }
}

std::string DEM::StoneMaterialContact::get_output_string() const {
    std::stringstream ss;
    ss << hmax_ << ", " << hs_ << ", " << hlinear_ << ", " << hp_ << ", " << FN_;
    return ss.str();
}

std::string DEM::StoneMaterialContact::restart_data() const {
    std::ostringstream ss;
    ss << named_print(FN_, "FN") << ", "
       << named_print(a_, "a") << ", "
       << named_print(R0_, "R0") << ", "
       << named_print(Fs_, "Fs") << ", "
       << named_print(hs_, "hs") << ", "
       << named_print(hlinear_, "hlinear") << ", "
       << named_print(h_, "h") << ", "
       << named_print(hmax_, "hmax") << ", "
       << named_print(Fmax_, "Fmax") << ", "
       << named_print(hp_, "hp") << ", "
       << named_print(hl_, "hl") << ", "
       << named_print(h1_, "h1") << ", "
       << named_print(FT_, "FT") << ", "
       << named_print(old_dT_, "old_dT") << ", "
       << named_print(load_var_tangential_, "load_var_tangential") << ", "
       << named_print(load_var_normal_, "load_var_normal") << ", "
       << named_print(turning_point_, "turning_point") << ", "
       << named_print(FT0[0], "FT0_1") << ", "
       << named_print(FT0[1], "FT0_2") << ", "
       << named_print(kp_, "kp") << ", "
       << named_print(ke_, "ke") << ", "
       << named_print(kl_, "kl") << ", "
       << named_print(ku_, "ku") << ", "
       << named_print(unloading_exp_, "unloading_exp") << ", "
       << named_print(mu_, "mu") << ", "
       << named_print(old_mu_, "old_mu") << ", "
       << named_print(kT_, "kT") << ", "
       << named_print(mu_r_, "mu_r") << ", "
       << named_print(kr_, "kr") << ", "
       << named_print(Tr_, "Tr") << ", "
       << named_print(fracture_strength_p1_, "fracture_strength_p1") << ", "
       << named_print(fracture_strength_p2_, "fracture_strength_p2") << ", "
       << named_print(id2_, "id2");
    return ss.str();
}
