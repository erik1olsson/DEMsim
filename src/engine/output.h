//
// Created by erolsson on 2018-07-26.
//

#ifndef DEMSIM_OUTPUT_H
#define DEMSIM_OUTPUT_H

#include <chrono>
#include <filesystem>
#include <string>
#include <vector>

#include "contact.h"
#include "../utilities/contact_matrix.h"
#include "../surfaces/surface_base.h"
#include "fine_cylinder.h"

namespace DEM {

    template<typename ForceModel, typename ParticleType>
    class Engine;
    class ParameterMap;

    template<typename ForceModel, typename ParticleType>
    class Output {
        using SurfaceType = Surface<ForceModel, ParticleType>;
        using ContactType = Contact<ForceModel, ParticleType>;;

    public:
        Output(const std::string& directory, std::chrono::duration<double> interval,
               const Engine<ForceModel, ParticleType>& engine, std::string name, bool remove_old_files=true);
        Output(const ParameterMap& parameters, const Engine<ForceModel, ParticleType>& engine);

        [[nodiscard]] const std::string& get_name() const { return name_; }

        void run_output();
        void run_output(const std::chrono::duration<double>& increment);

        bool print_particles = false;
        bool print_kinetic_energy = false;
        bool print_surface_positions = false;
        bool print_surface_forces = false;
        bool print_particle_cracks = false;
        bool print_contacts = false;
        bool print_bounding_box = false;
        bool print_periodic_bc = false;
        bool print_mirror_particles = false;
        bool print_fabric_force_tensor = false;
        bool print_fine_cylinder = false;
        bool print_fines_particle_forces = false;

        [[nodiscard]] std::string restart_data() const;
        void set_new_directory(const std::string& directory);
        void add_particle_to_follow(std::size_t particle_id);
        void activate() {active_ = true; }
        void deactivate() {active_ = false; }

    private:
        bool active_ = true;
        using OutputFunPtr = void (Output::*)() const;
        using FuncVec = std::vector<std::pair<bool&, OutputFunPtr>>;
        std::string name_;
        const std::vector<ParticleType*>& particles_;
        const std::vector<SurfaceType*>& surfaces_;
        const ContactMatrix<ContactType>& contacts_;
        const Engine<ForceModel, ParticleType>& engine_;

        const FineCylinder<ForceModel, ParticleType>* fine_cylinder_ = nullptr;

        std::filesystem::path directory_;
        std::vector<const ParticleType*> particles_to_print_;

        FuncVec output_functions_ {{print_particles,              &Output::write_particles},
                                   {print_kinetic_energy,         &Output::write_kinetic_energy},
                                   {print_surface_positions,      &Output::write_surface_positions},
                                   {print_surface_forces,         &Output::write_surface_forces},
                                   {print_particle_cracks,        &Output::write_particle_cracks},
                                   {print_contacts,               &Output::write_contacts},
                                   {print_bounding_box,           &Output::write_bounding_box},
                                   {print_periodic_bc,            &Output::write_periodic_bc},
                                   {print_mirror_particles,       &Output::write_mirror_particles},
                                   {print_fabric_force_tensor,    &Output::write_fabric_force_tensor},
                                   {print_fine_cylinder,          &Output::write_fine_cylinder},
                                   {print_fines_particle_forces,  &Output::write_fines_particle_forces}
        };

        std::chrono::duration<double> current_time_;
        std::chrono::duration<double> time_until_output_;
        std::chrono::duration<double> interval_;

        void write_particles() const;
        void write_kinetic_energy() const;
        void write_surface_positions() const;
        void write_surface_forces() const;
        void write_particle_cracks() const;
        void write_contacts() const;
        void write_bounding_box() const;
        void write_periodic_bc() const;
        void write_mirror_particles() const;
        void write_fabric_force_tensor() const;
        void write_particles_to_follow() const;
        void write_fine_cylinder() const;
        void write_fines_particle_forces() const;
    };
}

#include "output.tpp"

#endif //DEMSIM_OUTPUT_H
