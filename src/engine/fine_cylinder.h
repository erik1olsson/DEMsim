//
// Created by erolsson on 3/29/26.
//

#ifndef DEMSIM_FINECYLINDER_H
#define DEMSIM_FINECYLINDER_H

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

    private:
        FineMaterial* material_;
        double stone_density_;

        double p_ = 0;
        double D_;
        double Dideal_;
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
            double boundary_ratio = boundary_area/(3.1415*height()*radius_*radius_ - particle_volume);
            material_.update(D, boundary_ratio);
        }

        [[nodiscard]] double density() const {
            double slice_volume = 3.1415*height()*radius_*radius_;
            return fine_mass/(slice_volume - particle_volume);
        }
        [[nodiscard]] double relative_density() const {return material_.density();}
        [[nodiscard]] double ideal_density() const {return material_.ideal_density();}
        [[nodiscard]] double pressure() const {return material_.pressure(); }
        [[nodiscard]] double mass() const {return  fine_mass; }
        [[nodiscard]] double height() const {return  p2_.position - p1_.position; }
        [[nodiscard]] double radius() const {return  radius_; }
        [[nodiscard]] std::pair<double, double> get_positions() const { return {p1_.position, p2_.position}; }
        double particle_volume = 0;
        double boundary_area = 0;


        const FinesIntegrationPoint<ForceModel, ParticleType>& p1_;
        const FinesIntegrationPoint<ForceModel, ParticleType>& p2_;
    private:
        double fine_mass;
        double radius_;
        FineMaterialModel material_;

    };


    template<typename ParticleType>
    class ParticleForceModel {
    public:
        explicit ParticleForceModel(const ParticleType* p);
        [[nodiscard]] double get_force() const {return force_; }
        [[nodiscard]] const Vec3& get_displacement() const {return displacement_; }
        [[nodiscard]] double get_local_density() const {return local_density_; }
        void update(const Vec3& displacement, double force, double local_density) {
            displacement_ = displacement;
            force_ = force;
            local_density_ = local_density;
        }

        const ParticleType* get_particle() const { return particle_; }

    private:
        const ParticleType* particle_;
        Vec3 displacement_ = {0, 0,0};
        double force_ = 0;
        double local_density_;
    };

    template<typename ForceModel, typename ParticleType>
    class FineCylinder {
    public:
        FineCylinder(double radius, double height, double fine_mass, double stone_density, std::size_t numpoints,
            const std::vector<ParticleType*>& particles, FineMaterial* mat);
        void update(std::chrono::duration<double> dt, const Vec3& gravity);
        void connect_surface(PointSurface<ForceModel, ParticleType>* surface, std::size_t point);
        [[nodiscard]] std::string get_output_string() const;
        [[nodiscard]] double get_height() const;
        [[nodiscard]] double get_position_of_point(std::size_t) const;
        [[nodiscard]] double get_density_of_slice(std::size_t) const;
        [[nodiscard]] Vec3 get_particle_force(std::size_t id) const;


    private:
        double radius_;
        std::vector<FinesMaterialPoint<ForceModel, ParticleType>> mass_points_;
        std::vector<FinesIntegrationPoint<ForceModel, ParticleType>> integration_points_;
        std::vector<double> velocities_;
        std::vector<double> accelerations_;
        const std::vector<ParticleType*>& particles_;

        std::vector<ParticleForceModel<ParticleType>> particle_forces_;
        std::vector<std::pair<std::size_t, std::size_t> > particle_indices_;

        FineMaterial* material_;

    };
}

#include "fine_cylinder.tpp"

#endif //DEMSIM_FINECYLINDER_H