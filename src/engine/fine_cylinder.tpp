//
// Created by erolsson on 3/29/26.
//

#include <cmath>
#include <sstream>

#include "../surfaces/surface_base.h"

inline DEM::FineMaterial::FineMaterial(const double stone_density, const double initial_fine_density) :
    stone_density_(stone_density),
    D_(initial_fine_density/stone_density)
{

}

inline void DEM::FineMaterial::update(double density) {
    density /= stone_density_;
    double delta_d = density - D_;
    D_ = density;

    if (D_ > Dmax_ && D_ > D0) {
        // "Plastic" loading
        Dmax_ = D_;

        p_ = kl*pow(D_ - D0, n);
        pmax_ = p_;
        D0u = D0 + a*(Dmax_ - D0);
        ku = p_/pow(D_ - D0u, n);
    }

    else if (delta_d <= 0 && pmax_ > 0 && D_ > D0u) {
        // Unloading
        p_ = ku*pow(D_ - D0u, n);
        double gamma = pow(p_/pmax_, 1/(b*n));
        if (gamma < 1) {
            D0r = (D_ - gamma*Dmax_)/(1-gamma);
        }
        else {
            D0r = D0u;
        }
        kr = pmax_/pow(Dmax_ - D0r, b*n);
    }
    else if (delta_d > 0 && pmax_ > 0 && D_ > D0r) {
        // Reloading
        p_ = kr*pow(D_ - D0r, b*n);
    }
    else {
        p_ = 0;
    }
}

template<typename ForceModel, typename ParticleType>
DEM::FineCylinder<ForceModel, ParticleType>::FineCylinder(const double radius, double height, double fine_mass, double
                                                          stone_density, std::size_t numpoints)
    :radius_(radius){
    integration_points_.reserve(numpoints);
    integration_points_.emplace_back(0);
    velocities_ = std::vector<double>(numpoints, 0);
    accelerations_ = std::vector(numpoints, -9.82);
    const auto n_slices = static_cast<double>(numpoints -1 );
    for (std::size_t i = 1; i != numpoints; ++i) {
        double ip_pos = static_cast<double>(i)*height/n_slices;
        integration_points_.emplace_back(ip_pos);
        auto p = FinesMaterialPoint(integration_points_[i-1], integration_points_[i], fine_mass/n_slices,
            stone_density, radius);
        mass_points_.push_back(p);
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::FineCylinder<ForceModel, ParticleType>::update(std::chrono::duration<double> dt, const Vec3& gravity) {
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
        auto& int_point = integration_points_[i];
        if (int_point.surface != nullptr) {
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
           << mass_points_[i].pressure() << ", ";
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
