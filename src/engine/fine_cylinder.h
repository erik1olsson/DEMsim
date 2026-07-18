//
// Created by erolsson on 3/29/26.
//

#ifndef DEMSIM_FINECYLINDER_H
#define DEMSIM_FINECYLINDER_H

#include <chrono>
#include <string>
#include <vector>

#include "../surfaces/point_surface.h"

#include "../utilities/vec3.h"

class Particle;
class Cylinder;

namespace DEM {
    class FineMaterial {
    public:
        FineMaterial(double stone_density, double initial_fine_density);

        [[nodiscard]] double pressure() const {return p_; };
        [[nodiscard]] double density() const {return D_; };
        void update(double density);
    private:
        // Material parameters, hard coded for now
        constexpr static double D0 = 4.59801261e-01;
        constexpr static double kl = 3.16605936e+06*1e6;
        constexpr static double n = 4.20438353e+00;

        constexpr static double a = 0.658433977547588;
        constexpr static double b = 0.8;

        double stone_density_;

        double p_ = 0;
        double D_;

        // Internal state vars
        double Dmax_ = 0;

        // Convenience variables
        double pmax_ = 0;

        double D0u = D0;
        double ku = kl;

        double D0r = D0;
        double kr = kl;
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
            const double stone_density, const double radius) :
        p1_(p1), p2_(p2), fine_mass(fine_mass), radius_(radius), material_(stone_density, density())
        {}

        void update() {
            const double D = density();
            material_.update(D);
        }

        [[nodiscard]] double density() const { return fine_mass/(p2_.position - p1_.position)/3.1415/radius_/radius_;  }
        [[nodiscard]] double relative_density() const {return material_.density();}
        [[nodiscard]] double pressure() const {return material_.pressure(); }
        [[nodiscard]] double mass() const {return  fine_mass; }

    private:
        const FinesIntegrationPoint<ForceModel, ParticleType>& p1_;
        const FinesIntegrationPoint<ForceModel, ParticleType>& p2_;

        double fine_mass;
        double radius_;
        FineMaterial material_;
    };

    template<typename ForceModel, typename ParticleType>
    class FineCylinder {
    public:
        FineCylinder(double radius, double height, double fine_mass, double stone_density, std::size_t numpoints);
        void update(std::chrono::duration<double> dt, const Vec3& gravity);
        void connect_surface(PointSurface<ForceModel, ParticleType>* surface, std::size_t point);
        [[nodiscard]] std::string get_output_string() const;
        [[nodiscard]] double get_height() const;
        [[nodiscard]] double get_position_of_point(std::size_t) const;
        [[nodiscard]] double get_density_of_slice(std::size_t) const;


    private:
        double radius_;
        std::vector<FinesMaterialPoint<ForceModel, ParticleType>> mass_points_;
        std::vector<FinesIntegrationPoint<ForceModel, ParticleType>> integration_points_;
        std::vector<double> velocities_;
        std::vector<double> accelerations_;
    };
}

#include "fine_cylinder.tpp"

#endif //DEMSIM_FINECYLINDER_H