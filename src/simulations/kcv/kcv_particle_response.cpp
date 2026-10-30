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

void DEM::kcv_particle_response(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = StoneMaterialContact;
    using ParticleType = FractureableSphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;

    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto particle_file = parameters.get_parameter<std::string>("radius_file");
    auto gas_density = parameters.get_parameter<double>("gas_density");

    EngineType simulator(1us);

    auto mat = simulator.create_material<StoneMaterial>(2650.);
    mat->E = parameters.get_parameter<double>("E");
    mat->nu = parameters.get_parameter<double>("nu");
    mat->unloading_exponent = parameters.get_parameter<double>("unloading_exponent");
    mat->mu = parameters.get_parameter<double>("mu");
    mat->mu_wall = parameters.get_parameter<double>("mu_wall");

    auto compaction_velocity = parameters.get_parameter<double>("compaction_velocity");
    auto particle_radii = read_vector_from_file<double>(particle_file);

    double stone_volume = 0;
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

    simulator.create_cylinder(main_cylinder_radius, Vec3(0, 0, 1), Vec3(0, 0,0),
                              filling_height, true, true, false);

    // Creating The bottom plate surface
    Vec3 p1(-main_cylinder_radius, -main_cylinder_radius, 0.);
    Vec3 p2(-main_cylinder_radius,  main_cylinder_radius, 0.);
    Vec3 p3(main_cylinder_radius,   main_cylinder_radius, 0.);
    Vec3 p4(main_cylinder_radius,  -main_cylinder_radius, 0.);
    std::vector<Vec3> bottom_points{p4, p3, p2, p1};
    auto bottom_surface = simulator.create_point_surface(bottom_points, true);
    std::cout << "Normal of bottom surface is " << bottom_surface->get_normal() << "\n";

    // Creating The top plate surface
    Vec3 p5(-main_cylinder_radius, -main_cylinder_radius, filling_height);
    Vec3 p6(-main_cylinder_radius,  main_cylinder_radius, filling_height);
    Vec3 p7(main_cylinder_radius,   main_cylinder_radius, filling_height);
    Vec3 p8(main_cylinder_radius,  -main_cylinder_radius, filling_height);
    std::vector<Vec3> top_points{p5, p6, p7, p8};
    auto top_surface = simulator.create_point_surface(top_points, true);
    std::cout << "Normal of top surface is " << top_surface->get_normal() << "\n";

    std::cout << "Starting height is " << filling_height << "\n";

    auto particle_output = simulator.create_output(output_directory + "/animation", 0.01s);
    particle_output->print_particles = true;
    particle_output->print_surface_positions = true;
    particle_output->print_contacts = true;

    simulator.set_gravity(Vec3(0, 0, -9.820));
    simulator.set_mass_scale_factor(1.);
    std::vector<ParticleType*> particles;
    simulator.set_rotation(false);

    std::vector<double> simulation_particle_radii;
    simulation_particle_radii.assign(particle_radii.begin(), particle_radii.end());

    std::sort(simulation_particle_radii.begin(), simulation_particle_radii.end(), std::greater<>());

    std::cout << "Particles in this layer: " << particle_radii.end() -  particle_radii.begin() << "\n";

    auto particle_positions = random_fill_cylinder(0, filling_height,
                                                   main_cylinder_radius, simulation_particle_radii);

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

    auto output_shake = simulator.create_output(output_directory + "/shake/", 0.001s);
    output_shake->print_kinetic_energy = true;
    output_shake->print_surface_forces = true;
    output_shake->print_surface_positions = true;

    // Moving the lid to the top particle
    auto bbox = simulator.get_bounding_box();
    auto zmax1 = bbox[5];
    top_surface->move(Vec3(0, 0, zmax1 - filling_height), Vec3(0, 0, 0));
    std::cout << "Flipping gravity to up\n";
    simulator.set_gravity({0, 0, 9.82});
    run_for_time01.reset(0.1s);
    simulator.run(run_for_time01);
    max_velocity.set_new_value(0.1);
    simulator.run(max_velocity);

    // Moving the lid to the bottom particle
    bbox = simulator.get_bounding_box();
    auto zmin = bbox[4];
    bottom_surface->move(Vec3(0, 0, zmin), Vec3(0, 0, 0));
    std::cout << "Flipping gravity to down\n";
    simulator.set_gravity({0, 0, -9.82});
    run_for_time01.reset(0.1s);
    simulator.run(run_for_time01);
    max_velocity.set_new_value(0.1);
    simulator.run(max_velocity);

    // Dropping the large particle
    bbox = simulator.get_bounding_box();
    top_surface->move(Vec3(0, 0, filling_height - zmax1), Vec3(0, 0, 0));

    auto zmax2 = bbox[5];
    double large_radius = 0.008;
    simulator.create_particle(large_radius, Vec3(0, 0, zmax2 + 3*large_radius), Vec3(0, 0,0),
        mat);
    simulator.setup();
    output_shake->deactivate();
    auto output_drop = simulator.create_output(output_directory + "/particle_drop/", 0.001s);
    output_drop->print_kinetic_energy = true;
    output_drop->print_surface_forces = true;
    output_drop->print_surface_positions = true;

    // Moving the lid to the top particle
    run_for_time01.reset(0.1s);
    simulator.run(run_for_time01);
    max_velocity.set_new_value(0.1);
    simulator.run(max_velocity);

    // Compacting the particle down
    bbox = simulator.get_bounding_box();
    auto zmax3 = bbox[5];
    top_surface->move(Vec3(0, 0, zmax3 - filling_height), Vec3(0, 0, 0));
    top_surface->set_velocity(Vec3(0 ,0, -compaction_velocity));

    output_drop->deactivate();
    auto output_compaction = simulator.create_output(output_directory + "/compaction/", 0.001s);
    output_compaction->print_kinetic_energy = true;
    output_compaction->print_surface_forces = true;
    output_compaction->print_surface_positions = true;

    EngineType::RunForTime run_for_time5(simulator, 5s);
    simulator.run(run_for_time5);
}