//
// Created by erolsson on 12/08/2019.
//

#include <sstream>

#include "../simulations.h"

#include "../../engine/engine.h"

#include "../../contact_models/stone_material_contact.h"
#include "../../materials/stone_material.h"

#include "../../utilities/file_reading_functions.h"
#include "../../utilities/filling_functions.h"

#include "../../engine/fine_cylinder.h"

void DEM::kcv_particle_resistance(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = StoneMaterialContact;
    using ParticleType = FractureableSphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;

    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    const auto particle_file = parameters.get_parameter<std::string>("radius_file");
    const auto gas_density = parameters.get_parameter<double>("gas_density");

    EngineType simulator(1us);

    const auto mat = simulator.create_material<StoneMaterial>(2650.);
    mat->E = parameters.get_parameter<double>("E");
    mat->nu = parameters.get_parameter<double>("nu");
    mat->unloading_exponent = parameters.get_parameter<double>("unloading_exponent");
    mat->mu = parameters.get_parameter<double>("mu");
    mat->mu_wall = parameters.get_parameter<double>("mu_wall");

    auto compaction_velocity = parameters.get_parameter<double>("compaction_velocity");
    const auto particle_radii = read_vector_from_file<double>(particle_file);

    const auto test_particle_radius = parameters.get_parameter<double>("test_particle_radius");
    const auto test_particle_height = parameters.get_parameter<double>("test_particle_height");
    const auto pressure_levels = parameters.get_vector<double>("pressure_levels");
    const auto loading_time = parameters.get_parameter<double>("loading_time");
    const auto loading_distance = parameters.get_parameter<double>("loading_distance");
    auto particle_velocity = loading_distance/loading_time;
    double stone_volume = 4*3.1415*pow(test_particle_radius, 3)/3;
    for (const auto& r: particle_radii) {
        stone_volume += 4*3.1415*pow(r, 3)/3;
    }

    double stone_mass = stone_volume*mat->density;
    std::cout << "Mass of the stones: " << stone_mass << std::endl;
    double main_cylinder_radius = pow(stone_volume/3.1415, 1./3);
    double main_cylinder_height = 2*main_cylinder_radius;
    std::cout << "Target cylinder height: " << main_cylinder_height << std::endl;
    std::cout << "Target cylinder radius: " << main_cylinder_radius << std::endl;
    auto gas_volume = stone_volume/gas_density;

    auto filling_height = gas_volume/(3.1415*main_cylinder_radius*main_cylinder_radius);
    auto test_particle_z = test_particle_height*filling_height/main_cylinder_height;
    simulator.create_cylinder(main_cylinder_radius, Vec3(0, 0, 1), Vec3(0, 0,0),
                              filling_height, true, true, false);

    // Creating The bottom plate surface
    Vec3 p1(-main_cylinder_radius, -main_cylinder_radius, 0.);
    Vec3 p2(-main_cylinder_radius,  main_cylinder_radius, 0.);
    Vec3 p3(main_cylinder_radius,   main_cylinder_radius, 0.);
    Vec3 p4(main_cylinder_radius,  -main_cylinder_radius, 0.);
    std::vector bottom_points{p4, p3, p2, p1};
    auto bottom_surface = simulator.create_point_surface(bottom_points, true);
    std::cout << "Normal of bottom surface is " << bottom_surface->get_normal() << "\n";

    // Creating The top plate surface
    Vec3 p5(-main_cylinder_radius, -main_cylinder_radius, filling_height);
    Vec3 p6(-main_cylinder_radius,  main_cylinder_radius, filling_height);
    Vec3 p7(main_cylinder_radius,   main_cylinder_radius, filling_height);
    Vec3 p8(main_cylinder_radius,  -main_cylinder_radius, filling_height);
    std::vector top_points{p5, p6, p7, p8};
    auto top_surface = simulator.create_point_surface(top_points, true);
    std::cout << "Normal of top surface is " << top_surface->get_normal() << "\n";

    std::cout << "Starting height is " << filling_height << "\n";

    auto test_particle = simulator.create_particle(test_particle_radius, Vec3(0, 0, test_particle_z),
        Vec3(0,0,0), mat);
    std::vector<std::pair<Vec3, double>> placed_particles;
    placed_particles.emplace_back(test_particle->get_position(), test_particle->get_radius());

    auto particle_output = simulator.create_output(output_directory + "/animation", 0.01s);
    particle_output->print_particles = true;
    particle_output->print_surface_positions = true;

    simulator.set_gravity(Vec3(0, 0, -9.820));
    simulator.set_mass_scale_factor(1.);
    std::vector<ParticleType*> particles;
    simulator.set_rotation(false);

    std::vector<double> simulation_particle_radii;
    simulation_particle_radii.assign(particle_radii.begin(), particle_radii.end());

    std::sort(simulation_particle_radii.begin(), simulation_particle_radii.end(), std::greater<>());

    std::cout << "Particles in this layer: " << particle_radii.end() -  particle_radii.begin() << "\n";

    auto particle_positions = random_fill_cylinder(0, filling_height,
                                                   main_cylinder_radius, simulation_particle_radii, placed_particles);

    for (std::size_t i = 0; i != particle_positions.size(); ++i) {
        ParticleType *p = simulator.create_particle(simulation_particle_radii[i], particle_positions[i],
                                                    Vec3(0, 0, 0), mat);
        particles.push_back(p);
        double min_particle_mass = 4*3.1415*0.004*0.004*0.004/3*mat->density;
        if (p->get_mass() < min_particle_mass) {
            p->set_mass(min_particle_mass);
        }
    }

    auto output_filling = simulator.create_output(output_directory + "/filling/", 0.001s);
    output_filling->print_kinetic_energy = true;
    output_filling->print_surface_positions = true;
    output_filling->print_surface_forces = true;

    simulator.setup();
    EngineType::RunForTime run_for_time01(simulator, 0.1s);
    simulator.run(run_for_time01);


    EngineType::ParticleVelocityLess max_velocity(simulator, 0.1, 0.01s);
    simulator.run(max_velocity);
    simulator.remove_output(output_filling);
    // Moving the lid to the top particle

    auto bbox = simulator.get_bounding_box();
    auto zmax = bbox[5];

    top_surface->move(Vec3(0, 0, zmax - filling_height), Vec3(0, 0, 0));
    std::string restart_filename = "after_filling.res";
    simulator.write_restart_file(output_directory + "/" + restart_filename);

    for (const auto& pressure: pressure_levels) {
        // Compacting to the right pressure level
        double comp_force = pressure*1e6*3.1415*main_cylinder_radius*main_cylinder_radius;
        std::cout << "Compacting to pressure level p=" << pressure << " MPa, corresponding to a force of "
                  << comp_force << " N" << std::endl;
        auto simulator1 = EngineType(output_directory + "/" + restart_filename);
        auto output_compaction = simulator1.create_output(output_directory + "/compaction/", 0.001s);
        output_compaction->print_kinetic_energy = true;
        output_compaction->print_surface_forces = true;
        output_compaction->print_surface_positions = true;

        output_compaction->activate();

        EngineType::SurfaceNormalForceGreater surface_normal_force_greater(top_surface, comp_force);
        top_surface = dynamic_cast<decltype(top_surface)>(simulator1.get_surface(top_surface->get_id()));
        top_surface->set_velocity(Vec3(0 ,0, -compaction_velocity));
        simulator1.run(surface_normal_force_greater);
        std::stringstream restart_filename_ss;
        restart_filename_ss << "restart_p=" << pressure << "MPa.res";
        restart_filename = restart_filename_ss.str();
        top_surface->set_velocity(Vec3(0 ,0, 0.));
        simulator1.write_restart_file(output_directory + "/" + restart_filename);

        // Making the particle resistance test
        output_compaction->deactivate();
        std::stringstream output_name;
        output_name << "particle_resistance_p=" << pressure << "MPa";
        auto output_particle = simulator1.create_output(output_directory + "/" + output_name.str(), 0.001s);
        output_particle->print_kinetic_energy  = true;
        output_particle->print_surface_forces = true;
        output_particle->print_surface_positions = true;
        output_particle->add_particle_to_follow(test_particle->get_id());

        simulator1.exclude_particle_from_simulation(test_particle, Vec3(0, 0, -particle_velocity));
        EngineType::RunForTime run_for_time(simulator1, std::chrono::duration<double>(loading_time));
        simulator1.run(run_for_time);
        test_particle = simulator1.get_particle(test_particle->get_id());
        test_particle->set_velocity(Vec3(0,0, particle_velocity));
        run_for_time.reset(std::chrono::duration<double>(loading_time));
        simulator1.run(run_for_time);

    }
}