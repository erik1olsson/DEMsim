//
// Created by erolsson on 9/17/26.
//

#include <sstream>

#include "../simulations.h"

#include "../../engine/engine.h"

#include "../../contact_models/stone_material_contact.h"
#include "../../materials/stone_material.h"
#include "../../materials/fine_material.h"

#include "../../utilities/file_reading_functions.h"
#include "../../utilities/filling_functions.h"

#include "../../engine/fine_cylinder.h"

void DEM::kcv_particle_drop(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = StoneMaterialContact;
    using ParticleType = FractureableSphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;

    using namespace std::chrono_literals;
    double particle_radius = 0.008;
    SimulationParameters parameters(settings_file_name);
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto particle_file = parameters.get_parameter<std::string>("stone_radius_file");
    auto fines_file = parameters.get_parameter<std::string>("fines_radius_file");
    auto gas_density = parameters.get_parameter<double>("gas_density");

    EngineType simulator(1us);

    auto mat = simulator.create_material<StoneMaterial>(2650.);
    mat->E = parameters.get_parameter<double>("E");
    mat->nu = parameters.get_parameter<double>("nu");
    mat->unloading_exponent = parameters.get_parameter<double>("unloading_exponent");
    mat->mu = parameters.get_parameter<double>("mu");
    mat->mu_wall = parameters.get_parameter<double>("mu_wall");

    auto fine_material = simulator.create_material<FineMaterial>(2650);
    fine_material->assign_parameters(fuller5_25_parameters);

    auto fines_radii = read_vector_from_file<double>(fines_file);

    double fines_volume = 0;
    for (const auto& r: fines_radii) {
        fines_volume += 4*3.1415*pow(r, 3)/3;
    }

    double stone_volume = 0;

    // stone_volume *= 0;
    std::cout << "Fines volume is: " << fines_volume << std::endl;
    std::cout << "Stone volume is: " << stone_volume << std::endl;

    double main_cylinder_radius = pow((fines_volume + stone_volume)/3.1415, 1./3);
    double main_cylinder_height = 2*main_cylinder_radius;
    std::cout << "Target cylinder height: " << main_cylinder_height << std::endl;
    std::cout << "Target cylinder radius: " << main_cylinder_radius << std::endl;

    auto gas_volume = (fines_volume + stone_volume)/gas_density;

    auto filling_height = gas_volume/(3.1415*main_cylinder_radius*main_cylinder_radius);

    simulator.create_cylinder(main_cylinder_radius, Vec3(0, 0, 1), Vec3(0, 0,0),
                              filling_height, true, true, false);

    simulator.set_gravity(Vec3(0, 0, -9.820));
    std::vector<ParticleType*> particles;
    simulator.set_rotation(false);

    auto fine_cylinder = simulator.create_fine_cylinder(main_cylinder_radius, filling_height,
        fines_volume*mat->density,  mat->density, 25, fine_material);

    // Creating The bottom plate surface
    Vec3 p1(-main_cylinder_radius, -main_cylinder_radius, 0.);
    const Vec3 p2(-main_cylinder_radius,  main_cylinder_radius, 0.);
    Vec3 p3(main_cylinder_radius,   main_cylinder_radius, 0.);
    Vec3 p4(main_cylinder_radius,  -main_cylinder_radius, 0.);
    std::vector bottom_points{p4, p3, p2, p1};
    const auto bottom_surface = simulator.create_point_surface(bottom_points, true);
    fine_cylinder->connect_surface(bottom_surface, 0);
    std::cout << "Normal of bottom surface is " << bottom_surface->get_normal() << "\n";

    // Creating The top plate surface
    std::cout << "Starting height is " << filling_height << "\n";

    auto output_filling = simulator.create_output(output_directory + "/filling/", 0.001s);
    output_filling->print_kinetic_energy = true;
    output_filling->print_surface_positions = true;
    output_filling->print_surface_forces = true;
    output_filling->print_fine_cylinder = true;

    auto particle_output = simulator.create_output(output_directory + "/animation", 0.01s);
    particle_output->print_particles = true;
    particle_output->print_surface_positions = true;
    particle_output->print_contacts = true;

    simulator.setup();
    EngineType::RunForTime run_for_time1(simulator, 0.5s);
    simulator.run(run_for_time1);
    simulator.remove_output(output_filling);

    auto fine_height = fine_cylinder->get_height();

    auto output_drop = simulator.create_output(output_directory + "/particle_drop/", 0.001s);
    output_drop->print_kinetic_energy = true;
    output_drop->print_surface_forces = true;
    output_drop->print_surface_positions = true;
    output_drop->print_fine_cylinder = true;
    output_drop->print_fines_particle_forces = true;

    auto p = simulator.create_particle(particle_radius,
        Vec3(0, 0,fine_height + 3*particle_radius), Vec3(0, 0, 0), mat);
    auto z = fine_height + 5*particle_radius;
    Vec3 p5(-main_cylinder_radius, -main_cylinder_radius, z);
    Vec3 p6(-main_cylinder_radius,  main_cylinder_radius, z);
    Vec3 p7(main_cylinder_radius,   main_cylinder_radius, z);
    Vec3 p8(main_cylinder_radius,  -main_cylinder_radius, z);
    std::vector top_points{p5, p6, p7, p8};
    const auto top_surface = simulator.create_point_surface(top_points, true);


    simulator.setup();
    EngineType::ParticleVelocityLess max_velocity(simulator, 0.001, 0.01s);
    simulator.run(max_velocity);
    output_drop->deactivate();
    top_surface->move(Vec3(0, 0, p->get_position().z() + p->get_radius() - z), Vec3(0, 0, 0));
    top_surface->set_velocity(Vec3(0, 0, -2*particle_radius));
    auto output_compaction = simulator.create_output(output_directory + "/particle_compaction/", 0.001s);
    output_compaction->print_kinetic_energy = true;
    output_compaction->print_surface_forces = true;
    output_compaction->print_surface_positions = true;
    output_compaction->print_fine_cylinder = true;
    output_compaction->print_fines_particle_forces = true;

    EngineType::RunForTime run_for_time2(simulator, 1.s);
    simulator.run(run_for_time2);

}