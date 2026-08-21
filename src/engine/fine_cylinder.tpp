//
// Created by erolsson on 3/29/26.
//
#pragma once

#include <cmath>
#include <sstream>

inline DEM::FineMaterialModel::FineMaterialModel(const double stone_density, const double initial_fine_density,
    FineMaterial* mat) :
    material_(mat),
    stone_density_(stone_density),
    D_(initial_fine_density/stone_density),
    Dideal_(initial_fine_density/stone_density)

{

}

inline void DEM::FineMaterialModel::update(double density, double boundary_ratio) {
    density /= stone_density_;

    D_ = density;
    double delta_d = density/(1 - boundary_ratio*material_->db) - Dideal_;
    Dideal_ = density/(1 - boundary_ratio*material_->db);
    double D0 = material_->D_th;


    if (Dideal_ >= Dmax_ && Dideal_ > D0) {
        // "Plastic" loading
        p_ = material_->kl*pow(Dideal_ - D0, material_->n);
        pmax_ = p_;
        Dmax_ = Dideal_;
        D0u = D0 + material_->a*(Dmax_ - D0);
        ku = p_/pow(Dideal_ - D0u, material_->n);
    }

    else if (delta_d <= 0 && pmax_ > 0 && Dideal_ > D0u) {
        // Unloading
        p_ = ku*pow(Dideal_ - D0u, material_->n);
        double gamma = pow(p_/pmax_, 1/(material_->b*material_->n));
        if (gamma < 1) {
            D0r = (Dideal_ - gamma*Dmax_)/(1-gamma);
        }
        else {
            D0r = D0u;
        }
        kr = pmax_/pow(Dmax_ - D0r, material_->b*material_->n);
    }
    else if (delta_d > 0 && pmax_ > 0 && Dideal_ > D0r) {
        // Reloading
        p_ = kr*pow(Dideal_ - D0r, material_->b*material_->n);
    }
    else {
        p_ = 0;
    }
}

template<typename ParticleType>
DEM::ParticleForceModel<ParticleType>::ParticleForceModel(const ParticleType* p) :
particle_(p) {

}


template<typename ForceModel, typename ParticleType>
DEM::FineCylinder<ForceModel, ParticleType>::FineCylinder(const double radius, double height, double fine_mass, double
                                                          stone_density, std::size_t numpoints,
                                                          const std::vector<ParticleType*>& particles, FineMaterial* mat)
    :radius_(radius), particles_(particles), material_(mat){
    integration_points_.reserve(numpoints);
    integration_points_.emplace_back(0);
    velocities_ = std::vector<double>(numpoints, 0);
    accelerations_ = std::vector(numpoints, -9.82);;
    const auto n_slices = static_cast<double>(numpoints -1 );
    for (std::size_t i = 1; i != numpoints; ++i) {
        double ip_pos = static_cast<double>(i)*height/n_slices;
        integration_points_.emplace_back(ip_pos);
        auto p = FinesMaterialPoint(integration_points_[i-1], integration_points_[i], fine_mass/n_slices,
            stone_density, radius, mat);
        mass_points_.push_back(p);
    }
    for (const auto& p: particles_) {
        particle_indices_.emplace_back(p->get_id(), particle_forces_.size());
        particle_forces_.emplace_back(p);
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::FineCylinder<ForceModel, ParticleType>::update(std::chrono::duration<double> dt, const Vec3& gravity) {
    for (auto& mp: mass_points_) {
        mp.particle_volume = 0;
        mp.boundary_area = 2*3.1415*mp.radius()*mp.height();
        if (&mp == &mass_points_.front() || &mp == &mass_points_.back()) {
             mp.boundary_area += 3.1415*mp.radius()*mp.radius();
        }
    }

    // Distributing particle volumes and surface areas into fines slices
    for (const auto& p: particles_) {
        auto lowest_point = p->get_position().z() - p->get_radius();
        auto upper_point = p->get_position().z() + p->get_radius();
        auto mp_it = std::lower_bound(mass_points_.begin(), mass_points_.end(), lowest_point,
          [](const auto& mp, double value)
             {
                   return mp.get_positions().first < value;
             }
        );
        if (mp_it != mass_points_.begin()) {
            --mp_it;
        }

        double distributed_volume = 0;
        double distributed_area = 0;
        bool last_slice = false;
        while (!last_slice) {
            if ((mp_it + 1) == mass_points_.end() || (mp_it + 1)->get_positions().second > upper_point) {
                last_slice = true;
            }
            double h = mp_it->get_positions().second - lowest_point;
            double volume_below_slice = 3.1415*h*h/3*(3*p->get_radius() - h);
            double area_below_slice = 2*3.1415*p->get_radius()*h;
            mp_it->particle_volume += (volume_below_slice - distributed_volume);
            mp_it->boundary_area += (area_below_slice - distributed_area);

            distributed_volume += (volume_below_slice - distributed_volume);
            distributed_area += (area_below_slice - distributed_area);

            ++mp_it;
        }
        auto p_volume = 4*3.1415*pow(p->get_radius(), 3)/3;
        auto p_area = 4*3.1415*pow(p->get_radius(), 2);
        if (mp_it == mass_points_.end()) {
            --mp_it;
        }
        mp_it->particle_volume += (p_volume - distributed_volume);
        mp_it->boundary_area += (p_area - distributed_area);
    }

    // Updating particle forces from fines

    for (auto& p_force: particle_forces_) {
        auto particle_displacement = p_force.get_displacement();
        particle_displacement += p_force.get_particle()->get_displacement_this_increment();
        double du = particle_displacement.length() - p_force.get_displacement().length();
        auto movement_dir = particle_displacement.normal();
        auto fines_point = p_force.get_particle()->get_position() + p_force.get_particle()->get_radius()*movement_dir;
        auto ip_it = std::lower_bound(integration_points_.begin(), integration_points_.end(), fines_point.z(),
          [](const FinesIntegrationPoint<ForceModel, ParticleType>& ip, double value)
             {
                   return ip.position < value;
             }
        );
        if (ip_it != integration_points_.begin()) {
            ip_it -= 1;
        }

        auto active_slice = mass_points_[ip_it  - integration_points_.begin()];
        auto ideal_density = active_slice.ideal_density();
        auto particle_force = p_force.get_force();
        double local_density = p_force.get_local_density();
        if (p_force.get_local_density() == 0.) {
            local_density = ideal_density;
        }
        local_density += 30*du;

        if (local_density > material_->D_th && ideal_density > material_->D_th && local_density > ideal_density) {
            double dp = 0;
            if (local_density > ideal_density) {
                dp = material_->n*material_->kl*pow(local_density - material_->D_th, material_->n - 1)*30*du;
                double r = p_force.get_particle()->get_radius();
                particle_force += dp*3.1415*r*r;
                if (particle_force < 0) {
                    particle_force = 0;
                }
            }
            /*
            std::cout << "Computing particle force for particle" << p_force.get_particle()->get_id() << "\n";
            std::cout << "\tFines ideal density is: " << ideal_density << std::endl;
            std::cout << "\tLocal ideal density is: " << local_density << std::endl;
            std::cout << "\tdu: " << du << std::endl;
            std::cout << "\tpressure increase: " << dp << std::endl;
            std::cout << "\tMovement dir: " << movement_dir << std::endl;
            std::cout << "\tFines force on particle is: " << particle_force << std::endl;
            */
        }
        else {
            particle_force = 0;
            local_density = ideal_density;
            particle_displacement.set_zero();
        }
        p_force.update(particle_displacement, particle_force, local_density);
    }

    // Moving the slices
    for (auto& ip: integration_points_) {
        ip.force = 0;
        ip.mass = 0;
    }

    for (std::size_t i = 0; i!= mass_points_.size(); ++i) {
        mass_points_[i].update();
        integration_points_[i].force -= mass_points_[i].pressure()*3.1415*radius_*radius_;
        integration_points_[i].mass += mass_points_[i].mass()/2;
        integration_points_[i+1].force += mass_points_[i].pressure()*3.1415*radius_*radius_;
        integration_points_[i+1].mass += mass_points_[i].mass()/2;
    }
    for (std::size_t i = 0; i!= integration_points_.size(); ++i) {
        if (auto& int_point = integration_points_[i]; int_point.surface != nullptr) {
            int_point.position = int_point.surface->get_points()[0].z();
        }
        else {
            accelerations_[i] = gravity.z() + int_point.force/int_point.mass;
            velocities_[i] += accelerations_[i]*dt.count();
            int_point.position += (velocities_[i])*dt.count();

        }
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::FineCylinder<ForceModel, ParticleType>::connect_surface(PointSurface<ForceModel, ParticleType>* surface,
                                                                  std::size_t point) {
    auto& int_point = integration_points_[point];
    surface->set_fines_integration_point(&int_point);
    int_point.surface = surface;
}

template<typename ForceModel, typename ParticleType>
std::string DEM::FineCylinder<ForceModel, ParticleType>::get_output_string() const {
    std::ostringstream ss;
    for (std::size_t i = 0; i != mass_points_.size(); ++i) {
        ss << integration_points_[i].position << ", " << mass_points_[i].relative_density() << ", "
           << mass_points_[i].pressure() << ", " << mass_points_[i].particle_volume << ", "
           << mass_points_[i].boundary_area << ", ";
    }
    ss << integration_points_.back().position;
    return ss.str();
}

template<typename ForceModel, typename ParticleType>
double DEM::FineCylinder<ForceModel, ParticleType>::get_height() const {
    return integration_points_.back().position;
}

template<typename ForceModel, typename ParticleType>
double DEM::FineCylinder<ForceModel, ParticleType>::get_position_of_point(std::size_t i) const {
    return integration_points_[i].position;
}

template<typename ForceModel, typename ParticleType>
double DEM::FineCylinder<ForceModel, ParticleType>::get_density_of_slice(std::size_t i) const {
    return mass_points_[i].relative_density();

}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::FineCylinder<ForceModel, ParticleType>::get_particle_force(std::size_t id) const {
    auto idx_iter = std::lower_bound(particle_indices_.begin(), particle_indices_.end(), id,
        [](const auto& pair1, auto val) {return pair1.first < val; });
    double force = particle_forces_[idx_iter->second].get_force();
    const Vec3& disp = particle_forces_[idx_iter->second].get_displacement();
    if (disp.is_zero()) {
        return Vec3{0, 0, 0};
    }
    return force*disp.normal();
}

