//
// Created by erolsson on 2018-07-27.
//

#ifndef DEMSIM_SURFACE_H
#define DEMSIM_SURFACE_H

#include <array>
#include <memory>
#include <vector>

#include "../utilities/amplitude.h"
#include "../utilities/contact_vector.h"
#include "../utilities/contact_matrix.h"
#include "../utilities/vec3.h"


namespace DEM {
    template<typename ForceModel, typename ParticleType>
    class Contact;

    template<typename ForceModel, typename ParticleType>
    class Engine;

    class ParameterMap;

        template<typename ForceModel, typename ParticleType>
    class FinesIntegrationPoint;

    template<typename ForceModel, typename ParticleType>
    class Surface {
        using ContactType = Contact<ForceModel, ParticleType>;
        using ContactPointerType = typename ContactMatrix<ContactType>::PointerType;
        using ForceAmpPtr = std::shared_ptr<Amplitude>;


    public:
        double max_velocity{0.};
        bool force_control{ false };
        explicit Surface(std::size_t id, std::size_t collision_id, std::string  name, bool adhesive=false);
        explicit Surface(const ParameterMap& parameters);
        virtual ~Surface() = default;
        [[nodiscard]] std::size_t get_id() const { return object_id_; }
        [[nodiscard]] std::size_t get_collision_id() const { return collision_id_; }
        [[nodiscard]] const std::string& get_name() const { return name_; }
        [[nodiscard]] virtual Vec3 get_normal(const Vec3& position) const = 0;

        // virtual double get_curvature_radius() const = 0;  //  ToDo Implement later
        [[nodiscard]] virtual double distance_to_point(const Vec3& point) const = 0;
        [[nodiscard]] virtual Vec3 vector_to_point(const Vec3& point) const = 0;
        [[nodiscard]] virtual Vec3 get_displacement_this_increment(const Vec3& position) const = 0;

        virtual void move(const Vec3& distance, const Vec3& velocity) = 0;
        virtual void rotate(const Vec3& position, const Vec3& rotation_vector) = 0;


        [[nodiscard]] double get_mass() const { return mass_; }
        void set_mass(const double mass) {mass_ = mass; }
        [[nodiscard]] virtual std::string get_output_string() const = 0;
        [[nodiscard]] virtual const std::array<double, 6>& get_bounding_box_values() const { return bbox_values_; };
        [[nodiscard]] const Vec3& get_velocity() const { return velocity_; }
        [[nodiscard]] const Vec3& get_acceleration() const { return acceleration_; }

        void set_velocity(const Vec3& v) { velocity_ = v; }
        void set_acceleration(const Vec3& a) { acceleration_ = a; }

        [[nodiscard]] Vec3 get_tangential_displacement_this_inc(const Vec3& point) const;
        void rest();
        void set_adhesive(bool val) { adhesive_ = val; }
        [[nodiscard]] bool adhesive() const { return adhesive_; }

        std::vector<ParticleType*> get_contacting_particles() const;
        [[nodiscard]] double get_normal_force() const;
        [[nodiscard]] Vec3 get_tangential_force() const;
        [[nodiscard]] Vec3 get_total_force() const;
        void add_contact(ContactPointerType contact, std::size_t index_of_other_object);
        void remove_contact(std::size_t index_of_other_object);

        [[nodiscard]] const std::array<ForceAmpPtr, 3>& get_applied_forces() const { return force_control_amplitudes_; }
        void set_force_amplitude(const ForceAmpPtr& amplitude, char direction);
        void remove_force_amplitude(char direction);

        void set_fines_integration_point(const FinesIntegrationPoint<ForceModel, ParticleType>* p) {fines_integration_point_= p;}

        [[nodiscard]] virtual std::string restart_data() const;
        [[nodiscard]] virtual std::string type() const = 0;

    protected:
        Vec3 velocity_{ Vec3(0, 0, 0) };
        Vec3 acceleration_{ Vec3(0, 0, 0) };
        Vec3 displacement_this_inc_{ Vec3(0, 0, 0) };

        // Might be updated to std::vector<std::pair<Vec3, Vec3> >
        // to allow for multiple rotations
        Vec3 rotation_this_inc_{ Vec3(0, 0, 0) };
        Vec3 rotation_point_{ Vec3(0, 0, 0) };
        std::array<double, 6> bbox_values_{0, 0, 0,
                                           0, 0, 0};

        virtual void update_bounding_box() = 0;

    private:
        std::size_t object_id_;
        std::size_t collision_id_;
        std::string name_;
        double mass_ { 0 };
        bool adhesive_ {false};
        std::array<ForceAmpPtr, 3> force_control_amplitudes_ = {nullptr, nullptr, nullptr};

        ContactVector<ContactPointerType> contacts_{ContactVector<ContactPointerType>()};
        const FinesIntegrationPoint<ForceModel, ParticleType>* fines_integration_point_ = nullptr;
    };
}

#include "surface_base.tpp"
#endif //DEMSIM_SURFACE_H
