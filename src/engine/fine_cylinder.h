//
// Created by erolsson on 3/29/26.
//

#ifndef DEMSIM_FINECYLINDER_H
#define DEMSIM_FINECYLINDER_H

#include <algorithm>
#include <array>
#include <chrono>
#include <string>
#include <vector>

#include "../surfaces/point_surface.h"
#include "../materials/fine_material.h"

#include "../utilities/vec3.h"

class Particle;
class Cylinder;

namespace DEM {
    class FineMaterialModel {
    public:
        FineMaterialModel(double stone_density, double initial_fine_density, FineMaterial* mat);

        [[nodiscard]] double pressure() const {return p_; };
        [[nodiscard]] double density() const {return D_; };
        [[nodiscard]] double ideal_density() const {return Dideal_; };

        void update(double density, double boundary_ratio);

        [[nodiscard]] double stiffness() const {return stiffness_; }
        [[nodiscard]] double density_change() const { return dD_; }
        [[nodiscard]] double pressure_change() const {return dp_; }
        [[nodiscard]] std::array<double, 6> get_state() const {return {Dmax_, pmax_, D0u, ku, D0r, kr};}
        void set_state(const std::array<double, 6>& new_state);
        void copy_state_from(const FineMaterialModel& other);
        std::string state;

    private:
        FineMaterial* material_;
        double stone_density_;
        double stiffness_ = 0;
        double p_ = 0;
        double D_;
        double Dideal_ = 0;
        double dD_ideal_ = 0;
        double dD_ = 0;
        double dp_ = 0;
        // Internal state vars
        double Dmax_ = 0;

        // Convenience variables
        double pmax_ = 0;

        double D0u = 0;
        double ku = 0;

        double D0r = 0;
        double kr = 0;

    };

    template<typename ForceModel, typename ParticleType>
    class FinesIntegrationPoint {
    public:
        explicit FinesIntegrationPoint(double pos) : position(pos) {}
        double position;
        double mass = 0;
        double force = 0;
        double velocity = 0;
        PointSurface<ForceModel, ParticleType>* surface = nullptr;
    };

    template<typename ForceModel, typename ParticleType>
    class FinesMaterialPoint {
    public:
        FinesMaterialPoint(const FinesIntegrationPoint<ForceModel, ParticleType>& p1,
            const FinesIntegrationPoint<ForceModel, ParticleType>& p2, const double fine_mass,
            const double stone_density, const double radius, FineMaterial* mat) :
        p1_(p1), p2_(p2), fine_mass(fine_mass), radius_(radius), material_(stone_density, density(), mat)
        {}

        void update() {
            const double D = density();
            boundary_ratio = boundary_area/(3.1415*height()*radius_*radius_ - particle_volume);
            material_.update(D, boundary_ratio);
        }

        [[nodiscard]] double density() const {
            double slice_volume = 3.1415*height()*radius_*radius_;
            return fine_mass/(slice_volume - particle_volume);
        }
        [[nodiscard]] double relative_density() const {return material_.density();}
        [[nodiscard]] double ideal_density() const {return material_.ideal_density();}
        [[nodiscard]] double density_change() const {return material_.density_change();}
        [[nodiscard]] double pressure_change() const {return material_.pressure_change();}
        [[nodiscard]] double pressure() const {return material_.pressure(); }
        [[nodiscard]] double mass() const {return  fine_mass; }
        [[nodiscard]] double height() const {return  p2_.position - p1_.position; }
        [[nodiscard]] double center() const {return 0.5*(p1_.position + p2_.position); }
        [[nodiscard]] double radius() const {return  radius_; }
        [[nodiscard]] double material_velocity(double z) const {
            return p1_.velocity + (p2_.velocity - p1_.velocity)/height()*(z - p1_.position);
        }
        [[nodiscard]] std::array<double, 6> get_material_state() const {return  material_.get_state(); }
        [[nodiscard]] const FineMaterialModel& get_material_model() const {return material_; }
        [[nodiscard]] std::pair<double, double> get_positions() const { return {p1_.position, p2_.position}; }
        double particle_volume = 0;
        double boundary_area = 0;
        double boundary_ratio = 0;


        const FinesIntegrationPoint<ForceModel, ParticleType>& p1_;
        const FinesIntegrationPoint<ForceModel, ParticleType>& p2_;
    private:
        double fine_mass;
        double radius_;
        FineMaterialModel material_;

    };


    template<typename ForceModel, typename ParticleType>
    class ParticleForceModel {
    public:
        explicit ParticleForceModel(const ParticleType* p);
        [[nodiscard]] double get_force() const {return force_; }
        [[nodiscard]] const Vec3& get_displacement() const {return displacement_; }
        [[nodiscard]] double get_resistance_pressure() const {return resistance_pressure_; }
        [[nodiscard]] double get_background_pressure() const {return background_pressure_; }
        [[nodiscard]] double get_local_density() const {return local_density_; }
        [[nodiscard]] double get_mobilized_pressure() const {return mobilized_pressure_; }

        // Elastic local density increase produced by the magnitude of the 3-D
        // particle-fines relative displacement.  displacement_ and
        // mobilized_pressure_ are kept constitutively consistent after every
        // update: the pressure obtained from displacement_ is exactly the
        // stored mobilized pressure. Relative motion that cannot be stored in
        // that state is plastic flow around the particle.
        [[nodiscard]] static constexpr double density_gradient() {
            return 45.0;
        }

        [[nodiscard]] static double density_increment(double displacement) {
            return density_gradient()*displacement;
        }

        [[nodiscard]] static double displacement_from_density_increment(
            double density_increment) {
            return std::max(0.0, density_increment)/density_gradient();
        }

        void update(const Vec3& displacement, double force, double local_density,
                    double resistance_pressure, double mobilized_pressure,
                    double background_pressure) {
            displacement_ = displacement;
            force_ = force;
            local_density_ = local_density;
            resistance_pressure_ = resistance_pressure;
            mobilized_pressure_ = mobilized_pressure;
            background_pressure_ = background_pressure;
        }

        const ParticleType* get_particle() const { return particle_; }
        FinesMaterialPoint<ForceModel, ParticleType>* active_slice = nullptr;

    private:
        const ParticleType* particle_;
        Vec3 displacement_ = {0, 0,0};
        double force_ = 0;
        double local_density_ = 0;
        double mobilized_pressure_ = 0;
        double resistance_pressure_ = 0;
        double background_pressure_ = 0;
    };

    template<typename ForceModel, typename ParticleType>
    class FineCylinder {
    public:
        FineCylinder(double radius, double height, double fine_mass, double stone_density, std::size_t numpoints,
            const std::vector<ParticleType*>& particles, FineMaterial* mat,
            double particle_yield_density_increment = 0.005);
        void add_particle(const ParticleType* particle);
        void update(std::chrono::duration<double> dt, const Vec3& gravity);
        void connect_surface(PointSurface<ForceModel, ParticleType>* surface, std::size_t point);
        [[nodiscard]] std::string get_output_string() const;
        [[nodiscard]] std::vector<std::string> particle_force_strings() const;
        [[nodiscard]] double get_height() const;
        [[nodiscard]] double get_kinetic_energy() const {return kinetic_energy_; }
        [[nodiscard]] double get_position_of_point(std::size_t) const;
        [[nodiscard]] double get_density_of_slice(std::size_t) const;
        [[nodiscard]] Vec3 get_particle_force(std::size_t id) const;


    private:
        struct LocalFinesState {
            double density = 0;
            double pressure = 0;
            double boundary_ratio = 0;
        };

        // The material variables are stored as slice averages.  For the
        // particle interaction they are interpreted as values at the slice
        // centres and interpolated linearly in between.  Outside the first and
        // last centre the closest value is kept constant; no extrapolation is
        // performed towards the cylinder boundaries.
        [[nodiscard]] LocalFinesState local_fines_state(double z) const;

        double radius_;
        std::vector<FinesMaterialPoint<ForceModel, ParticleType>> mass_points_;
        std::vector<FinesIntegrationPoint<ForceModel, ParticleType>> integration_points_;
        std::vector<double> velocities_;
        std::vector<double> accelerations_;
        const std::vector<ParticleType*>& particles_;

        std::vector<ParticleForceModel<ForceModel, ParticleType>> particle_forces_;
        std::vector<std::pair<std::size_t, std::size_t> > particle_indices_;

        FineMaterial* material_;

        // Maximum additional elastic density increment after the local fines
        // have reached the load-carrying density. This is the single new
        // particle-resistance parameter. Excess relative motion is plastic flow.
        double particle_yield_density_increment_;
        double kinetic_energy_ = 0;

    };
}

#include "fine_cylinder.tpp"

#endif //DEMSIM_FINECYLINDER_H
