//
// Created by erolsson on 3/29/26.
//
#pragma once

#include <cmath>
#include <iostream>
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
    dD_ = density - D_;
    D_ = density;

    double boundary_correction = (1 - boundary_ratio*material_->db);
    dD_ideal_ = density/boundary_correction - Dideal_;
    Dideal_ = density/boundary_correction;
    double D0 = material_->D_th;
    double new_pressure = 0;

    if (Dideal_ >= Dmax_ && Dideal_ > D0) {
        // "Plastic" loading
        new_pressure = material_->kl*pow((Dideal_ - D0)*boundary_correction, material_->n);
        stiffness_ = material_->n*material_->kl*pow((Dideal_ - D0)*boundary_correction, material_->n - 1);
        pmax_ = new_pressure;
        Dmax_ = Dideal_;
        D0u = D0 + material_->a*(Dmax_ - D0);
        ku = new_pressure/pow((Dideal_ - D0u)*boundary_correction, material_->n);
        state = "plastic loading";
    }

    else if (dD_ideal_ <= 0 && pmax_ > 0 && Dideal_ > D0u) {
        // Unloading
        state = "unloading";
        new_pressure = ku*pow((Dideal_ - D0u)*boundary_correction, material_->n);
        stiffness_ = material_->n*ku*pow((Dideal_ - D0u)*boundary_correction, material_->n - 1);
        double gamma = pow(new_pressure/pmax_, 1/(material_->b*material_->n));
        if (gamma < 1) {
            D0r = (Dideal_ - gamma*Dmax_)/(1-gamma);
        }
        else {
            D0r = D0u;
        }
        kr = pmax_/pow((Dmax_ - D0r)*boundary_correction, material_->b*material_->n);
    }
    else if (dD_ideal_ > 0 && pmax_ > 0 && Dideal_ > D0r) {
        // Reloading
        state = "reloading";
        new_pressure = kr*pow((Dideal_ - D0r)*boundary_correction, material_->b*material_->n);
        stiffness_ = material_->b*material_->n*kr*pow((Dideal_ - D0r)*boundary_correction, material_->b*material_->n - 1);
        ku = new_pressure/pow((Dideal_ - D0u)*boundary_correction, material_->n);
    }
    else {
        new_pressure = 0;
        stiffness_ = 0;
    }
    dp_ = new_pressure - p_;
    p_ = new_pressure;
}

inline void DEM::FineMaterialModel::set_state(const std::array<double, 6>& new_state) {
    Dmax_ = new_state[0];
    pmax_ = new_state[1];
    D0u = new_state[2];
    ku = new_state[3];
    D0r = new_state[4];
    kr = new_state[5];
}

inline void DEM::FineMaterialModel::copy_state_from(const FineMaterialModel& other) {
    stiffness_ = other.stiffness_;
    p_ = other.p_;
    D_ = other.D_;
    Dideal_ = other.Dideal_;
    dD_ideal_ = other.dD_ideal_;
    dD_ = other.dD_;
    dp_ = other.dp_;

    Dmax_ = other.Dmax_;
    pmax_ = other.pmax_;
    D0u = other.D0u;
    ku = other.ku;
    D0r = other.D0r;
    kr = other.kr;
    state = other.state;
}


template<typename ForceModel, typename ParticleType>
DEM::ParticleForceModel<ForceModel, ParticleType>::ParticleForceModel(const ParticleType* p) :
particle_(p) {

}


template<typename ForceModel, typename ParticleType>
DEM::FineCylinder<ForceModel, ParticleType>::FineCylinder(const double radius, double height, double fine_mass, double
                                                          stone_density, std::size_t numpoints,
                                                          const std::vector<ParticleType*>& particles, FineMaterial* mat,
                                                          const double particle_yield_density_increment)
    :radius_(radius), particles_(particles), material_(mat),
     particle_yield_density_increment_(std::max(0.0, particle_yield_density_increment)) {
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
       add_particle(p);
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::FineCylinder<ForceModel, ParticleType>::add_particle(const ParticleType* particle) {
    particle_indices_.emplace_back(particle->get_id(), particle_forces_.size());
    particle_forces_.emplace_back(particle);
}

template<typename ForceModel, typename ParticleType>
typename DEM::FineCylinder<ForceModel, ParticleType>::LocalFinesState
DEM::FineCylinder<ForceModel, ParticleType>::local_fines_state(const double z) const {
    const auto state_at = [this](const auto& slice) {
        return LocalFinesState{
            slice.density()/material_->density,
            slice.pressure(),
            slice.boundary_ratio
        };
    };

    if (mass_points_.empty()) {
        return {};
    }

    // Constant continuation from the end slice centres to the boundaries.
    // This avoids inventing a gradient outside the region supported by two
    // neighbouring slice values.
    if (mass_points_.size() == 1 || z <= mass_points_.front().center()) {
        return state_at(mass_points_.front());
    }
    if (z >= mass_points_.back().center()) {
        return state_at(mass_points_.back());
    }

    // First slice centre at or above z. The preceding slice is then the lower
    // interpolation point. Slice centres remain ordered as long as the 1-D
    // fines mesh has not inverted.
    const auto upper = std::lower_bound(
        mass_points_.begin() + 1, mass_points_.end(), z,
        [](const auto& slice, const double value) {
            return slice.center() < value;
        });
    const auto lower = upper - 1;

    const double lower_z = lower->center();
    const double upper_z = upper->center();
    const double centre_distance = upper_z - lower_z;
    if (!(centre_distance > 0)) {
        // A non-positive distance means that the continuum mesh has inverted.
        // Do not divide by zero here; the mesh dynamics can handle/report the
        // underlying problem separately.
        return state_at(*lower);
    }

    const double weight = std::clamp(
        (z - lower_z)/centre_distance, 0.0, 1.0);
    const auto lower_state = state_at(*lower);
    const auto upper_state = state_at(*upper);

    const auto interpolate = [weight](const double a, const double b) {
        return a + weight*(b - a);
    };

    return {
        interpolate(lower_state.density, upper_state.density),
        interpolate(lower_state.pressure, upper_state.pressure),
        interpolate(lower_state.boundary_ratio, upper_state.boundary_ratio)
    };
}

template<typename ForceModel, typename ParticleType>
void DEM::FineCylinder<ForceModel, ParticleType>::update(std::chrono::duration<double> dt, const Vec3& gravity) {
    for (auto& ip: integration_points_) {
        ip.force = 0;
        ip.mass = 0;
    }
    for (auto& mp: mass_points_) {
        mp.particle_volume = 0;
        mp.boundary_area = 2*3.1415*mp.radius()*mp.height();
        if ((&mp == &mass_points_.front() && mass_points_.front().p1_.surface ) ||
            (&mp == &mass_points_.back() && mass_points_.back().p2_.surface) ) {
             mp.boundary_area += 3.1415*mp.radius()*mp.radius();
        }
    }

    // Distributing particle volumes and surface areas into fines slices
    for (const auto& p: particles_) {
        auto lowest_point = p->get_position().z() - p->get_radius();
        if (lowest_point > get_height()) {
            continue;
        }
        auto upper_point = p->get_position().z() + p->get_radius();

        auto mp_it = std::lower_bound(mass_points_.begin(), mass_points_.end(), lowest_point,
          [](const auto& mp, double value)
             {
                   return mp.get_positions().first < value;
             }
        );

        if (lowest_point < 0) {
            mp_it = mass_points_.begin();
        }

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
            if (mp_it->get_positions().second > upper_point) {
                h = 2*p->get_radius();
            }
            double volume_below_slice = 3.1415*h*h/3*(3*p->get_radius() - h);
            double area_below_slice = 2*3.1415*p->get_radius()*h;
            mp_it->particle_volume += (volume_below_slice - distributed_volume);
            if (mp_it->particle_volume < 0) {
                std::abort();
            }
            mp_it->boundary_area += (area_below_slice - distributed_area);
            distributed_volume += (volume_below_slice - distributed_volume);
            distributed_area += (area_below_slice - distributed_area);
            ++mp_it;
        }
        auto p_volume = 4*3.1415*pow(p->get_radius(), 3)/3;
        auto p_area = 4*3.1415*pow(p->get_radius(), 2);
        if (mp_it == mass_points_.end()) {
            continue;
        }
        mp_it->particle_volume += (p_volume - distributed_volume);
        if (mp_it->particle_volume < 0) {
            std::abort();
        }
        mp_it->boundary_area += (p_area - distributed_area);
    }

    // Updating particle forces from fines

    for (std::size_t i = 0; i!= mass_points_.size(); ++i) {
        mass_points_[i].update();
    }

    for (auto& p_force: particle_forces_) {
        const Vec3 old_displacement = p_force.get_displacement();

        const double particle_radius = p_force.get_particle()->get_radius();
        const auto& particle_position = p_force.get_particle()->get_position();
        /*
        std::cout << "The mass points:\n";
        for (const auto& mp: mass_points_) {
            std::cout << "\tz1: " <<  mp.get_positions().first << ", z2: " << mp.get_positions().second
                      << ", pressure: " << mp.pressure() << " ideal density: " << mp.ideal_density()
                      << ", particle volume: " << mp.particle_volume << ", boundary area: " << mp.boundary_area << "\n";
        }
        */
        const double particle_lowest_point = particle_position.z() - particle_radius;
        const double particle_highest_point = particle_position.z() + particle_radius;
        const bool overlaps_fines = particle_highest_point > 0 &&
                                    particle_lowest_point < get_height();

        if (overlaps_fines) {
            // The 1-D fines state is sampled at the particle centre. Using the
            // leading particle surface would move the sampling point by one
            // diameter at every reversal and introduce a discontinuous change
            // of slice and material history.
            const double z_tolerance = 1.e-12*std::max(1.0, get_height());
            const double fines_z = std::clamp(
                particle_position.z(), z_tolerance, get_height() - z_tolerance);

            // Locate the slice every step. The compacted zone moves with the
            // particle and must not remain connected to a slice that it has
            // already left.
            auto slice_it = std::lower_bound(mass_points_.begin(), mass_points_.end(), fines_z,
              [](const FinesMaterialPoint<ForceModel, ParticleType>& mp, double value)
                 {
                       return mp.get_positions().second < value;
                 }
            );
            const auto active_slice_index = static_cast<std::size_t>(
                slice_it - mass_points_.begin());
            auto* new_active_slice = &mass_points_[active_slice_index];

            if (p_force.active_slice != new_active_slice) {
                p_force.active_slice = new_active_slice;
            }
            auto active_slice = p_force.active_slice;

            // Density, background pressure and boundary correction are smooth
            // local fields obtained from the values at the slice centres. The
            // active slice is used only for velocity interpolation and reaction
            // force distribution; no second hysteretic material history is
            // attached to the particle.
            const auto local_fines = local_fines_state(fines_z);
            const double slice_density = local_fines.density;
            const double background_pressure = local_fines.pressure;
            const double local_boundary_ratio = local_fines.boundary_ratio;

            // Split particle-fines relative motion into a bounded elastic part
            // and plastic flow around the particle. Loose fines may first be
            // compacted to the local load-carrying density. Beyond that point,
            // only particle_yield_density_increment_ is stored elastically.
            const double boundary_correction = std::max(
                1.e-12, 1.0 - local_boundary_ratio*material_->db);
            const double onset_density = material_->D_th*boundary_correction;
            const double density_increment_to_onset = std::max(
                0.0, onset_density - slice_density);
            const double maximum_density_increment =
                density_increment_to_onset + particle_yield_density_increment_;
            const double maximum_elastic_displacement =
                p_force.displacement_from_density_increment(
                    maximum_density_increment);

            // The continuum fields vary only with z, but the local particle
            // resistance is three-dimensional. Fines have only a vertical bulk
            // velocity; their unresolved horizontal velocity is zero in this
            // reduced model.
            const auto limit_elastic_displacement =
                [maximum_elastic_displacement](Vec3 displacement) {
                    const double length = displacement.length();
                    if (length > maximum_elastic_displacement && length > 0) {
                        displacement *= maximum_elastic_displacement/length;
                    }
                    return displacement;
                };
            const Vec3 bounded_old_displacement =
                limit_elastic_displacement(old_displacement);
            const double fines_velocity =
                active_slice->material_velocity(fines_z);
            Vec3 relative_velocity = p_force.get_particle()->get_velocity();
            relative_velocity -= Vec3(0, 0, fines_velocity);
            Vec3 trial_displacement = bounded_old_displacement;
            trial_displacement += relative_velocity*dt.count();
            Vec3 particle_displacement =
                limit_elastic_displacement(trial_displacement);
            const double trial_length = trial_displacement.length();
            const double plastic_flow_increment = std::max(
                0.0, trial_length -
                     maximum_elastic_displacement);

            const double old_length = bounded_old_displacement.length();
            const double new_length = particle_displacement.length();
            const double old_density_increment =
                p_force.density_increment(old_length);
            const double new_density_increment =
                p_force.density_increment(new_length);

            // Use only the continuous loading envelope for the particle-local
            // pressure increment. Hysteresis remains in the continuum material
            // model; using another complete FineMaterialModel here caused
            // discontinuous pressure jumps when its branch changed.
            const auto envelope_pressure =
                [this, onset_density](const double density) {
                    const double load_carrying_density = std::max(
                        0.0, density - onset_density);
                    return material_->kl*std::pow(
                        load_carrying_density, material_->n);
                };
            const double envelope_background_pressure =
                envelope_pressure(slice_density);
            const auto incremental_pressure =
                [&envelope_pressure, envelope_background_pressure,
                 slice_density](const double density_increment) {
                    return std::max(
                        0.0,
                        envelope_pressure(slice_density + density_increment) -
                        envelope_background_pressure);
                };

            const double old_available_pressure =
                incremental_pressure(old_density_increment);
            const double available_pressure =
                incremental_pressure(new_density_increment);

            // The unbounded trial state is used only to mobilize pressure
            // smoothly during continued plastic flow. Limit its evaluation to
            // one extra yield-density interval for numerical robustness; the
            // applied pressure is still capped by available_pressure below.
            const double trial_density_increment = std::min(
                p_force.density_increment(trial_length),
                maximum_density_increment + particle_yield_density_increment_);
            const double trial_available_pressure =
                incremental_pressure(trial_density_increment);

            // Background compaction may raise the available capacity but is
            // not allowed to increase the already mobilized pressure.
            // Background unloading may lower the force through the capacity
            // cap. Only relative particle-fines motion changes the mobilized
            // pressure.
            double mobilized_pressure = std::clamp(
                p_force.get_mobilized_pressure(),
                0.0, old_available_pressure);
            double movement_pressure_increment = 0;

            const double old_trial_dot =
                bounded_old_displacement.x()*trial_displacement.x() +
                bounded_old_displacement.y()*trial_displacement.y() +
                bounded_old_displacement.z()*trial_displacement.z();
            const bool direction_reversal =
                old_length > 1.e-15 && trial_length > 1.e-15 &&
                old_trial_dot < -0.99*old_length*trial_length;
            if (direction_reversal) {
                // A nearly collinear reversal has crossed the origin between
                // two time levels: unload the old compacted zone before the
                // opposite zone is mobilized. General 3-D rotations remain
                // continuous and do not trigger this branch.
                mobilized_pressure = 0;
                movement_pressure_increment = std::min(
                    trial_available_pressure, available_pressure);
            }
            else if (trial_length > old_length) {
                movement_pressure_increment = std::max(
                    0.0,
                    trial_available_pressure - old_available_pressure);
            }
            else if (new_length < old_length) {
                movement_pressure_increment = std::min(
                    0.0,
                    available_pressure - old_available_pressure);
            }

            mobilized_pressure = std::clamp(
                mobilized_pressure + movement_pressure_increment,
                0.0, available_pressure);
            double resistance_pressure = mobilized_pressure;
            const double new_local_pressure =
                background_pressure + available_pressure;

            double local_density = slice_density + new_density_increment;
            if (new_length <= 1.e-15) {
                // No stored relative deformation means no directional force.
                particle_displacement = Vec3(0, 0, 0);
                local_density = slice_density;
                resistance_pressure = 0;
                mobilized_pressure = 0;
            }

            const double r = p_force.get_particle()->get_radius();
            const double particle_force = 3.1415*r*r*resistance_pressure;
            p_force.update(particle_displacement, particle_force, local_density,
                           resistance_pressure, mobilized_pressure,
                           background_pressure);

            // The caller applies -particle_force*direction to the particle.
            // Apply the equal and opposite vertical reaction to the 1-D fines
            // continuum and distribute it to the two slice integration points
            // using the same linear shape functions as the velocity field. The
            // horizontal reaction is carried by the unresolved fines skeleton
            // and cylinder boundary; it is reported below but has no continuum
            // degree of freedom to which it can be applied.
            Vec3 fines_reaction_force(0, 0, 0);
            if (particle_force > 0 && !particle_displacement.is_zero()) {
                fines_reaction_force =
                    particle_force*particle_displacement.normal();

                const double slice_height = active_slice->height();
                if (slice_height > 0) {
                    const double upper_weight = std::clamp(
                        (fines_z - active_slice->p1_.position)/slice_height,
                        0.0, 1.0);
                    const double lower_weight = 1.0 - upper_weight;
                    integration_points_[active_slice_index].force +=
                        lower_weight*fines_reaction_force.z();
                    integration_points_[active_slice_index + 1].force +=
                        upper_weight*fines_reaction_force.z();
                }
            }

            if (false) {
                std::cout << p_force.get_particle()->get_id() << ", "
                          << p_force.get_particle()->get_position().z() << ", "
                          << particle_displacement.z() << ", "
                          << fines_z << ", "
                          << active_slice->p1_.position << ", "
                          << active_slice->p2_.position << ", "
                          << fines_velocity << ", "
                          << p_force.get_particle()->get_velocity().z() << ", "
                          << active_slice->p1_.velocity << ", "
                          << active_slice->p2_.velocity << ", "
                          << resistance_pressure << ", "
                          << mobilized_pressure << ", "
                          << available_pressure << ", "
                          << background_pressure << ", "
                          << new_local_pressure << ", "
                          << new_density_increment << ", "
                          << local_density << ", "
                          << slice_density << ", "
                          << new_density_increment - old_density_increment << ", "
                          << particle_force << ", "
                          << active_slice->density_change() << ", "
                          << maximum_density_increment << ", "
                          << plastic_flow_increment << ", "
                          << fines_reaction_force.z() << ", "
                          << envelope_background_pressure << ", "
                          << movement_pressure_increment << ", "
                          << particle_displacement.x() << ", "
                          << particle_displacement.y() << ", "
                          << relative_velocity.x() << ", "
                          << relative_velocity.y() << ", "
                          << fines_reaction_force.x() << ", "
                          << fines_reaction_force.y() << "\n";
            }
        }

        else {
            p_force.update(Vec3(0,0, 0), 0, 0, 0, 0, 0);
            p_force.active_slice = nullptr;
        }
    }

    // Moving the slices
    for (std::size_t i = 0; i!= mass_points_.size(); ++i) {
        integration_points_[i].force -= mass_points_[i].pressure()*3.1415*radius_*radius_;
        integration_points_[i].mass += mass_points_[i].mass()/2;
        integration_points_[i+1].force += mass_points_[i].pressure()*3.1415*radius_*radius_;
        integration_points_[i+1].mass += mass_points_[i].mass()/2;
    }
    kinetic_energy_ = 0;
    for (std::size_t i = 0; i!= integration_points_.size(); ++i) {
        if (auto& int_point = integration_points_[i]; int_point.surface != nullptr) {
            int_point.position = int_point.surface->get_points()[0].z();
        }
        else {
            accelerations_[i] = gravity.z() + int_point.force/int_point.mass;
            velocities_[i] += accelerations_[i]*dt.count();
            int_point.velocity = velocities_[i];
            kinetic_energy_ += int_point.mass*velocities_[i]*velocities_[i]/2;
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
std::vector<std::string> DEM::FineCylinder<ForceModel, ParticleType>::particle_force_strings() const {
    std::vector<std::string> strings;
    for (const auto& p_force: particle_forces_) {
        std::stringstream ss;
        const ParticleType* p = p_force.get_particle();
        ss << p->get_id() << ", " << p->get_position().x() << ", " << p->get_position().y() << ", "
           << p->get_position().z() << ", " << p_force.get_local_density() << ", " << p_force.get_force() << ", "
           << p_force.get_displacement().x() << ", " << p_force.get_displacement().y() << ", "
           << p_force.get_displacement().z();
        strings.push_back(ss.str());
    }
    return strings;
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
