//
// Created by erolsson on 12/08/2019.
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

void DEM::kcv_fines(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = StoneMaterialContact;
    using ParticleType = FractureableSphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;

    using namespace std::chrono_literals;

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

    auto compaction_velocity = parameters.get_parameter<double>("compaction_velocity");
    auto particle_radii = read_vector_from_file<double>(particle_file);
    auto fines_radii = read_vector_from_file<double>(fines_file);

    double fines_volume = 0;
    for (const auto& r: fines_radii) {
        fines_volume += 4*3.1415*pow(r, 3)/3;
    }

    double stone_volume = 0;
    for (const auto& r: particle_radii) {
        stone_volume += 4*3.1415*pow(r, 3)/3;
    }
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

    std::vector<double> simulation_particle_radii;
    simulation_particle_radii.assign(particle_radii.begin(), particle_radii.end());

    std::sort(simulation_particle_radii.begin(), simulation_particle_radii.end(), std::greater<>());

    std::cout << "Particles in this layer: " << particle_radii.end() -  particle_radii.begin() << "\n";

    auto particle_positions = random_fill_cylinder(0, filling_height,
                                                   main_cylinder_radius, simulation_particle_radii);
    //for (std::size_t i = 0; i != 1; ++i) {
    for (std::size_t i = 0; i != particle_positions.size(); ++i) {
        ParticleType *p = simulator.create_particle(simulation_particle_radii[i], particle_positions[i],
                                                    Vec3(0, 0, 0), mat);
        particles.push_back(p);
        if (double min_particle_mass = 4*3.1415*0.004*0.004*0.004/3*mat->density; p->get_mass() < min_particle_mass) {
            p->set_mass(min_particle_mass);
        }
    }

    auto fine_cylinder = simulator.create_fine_cylinder(main_cylinder_radius, filling_height,
        fines_volume*mat->density,  mat->density, 25, fine_material);


    // Creating The bottom plate surface
    Vec3 p1(-main_cylinder_radius, -main_cylinder_radius, 0.);
    const Vec3 p2(-main_cylinder_radius,  main_cylinder_radius, 0.);
    Vec3 p3(main_cylinder_radius,   main_cylinder_radius, 0.);
    Vec3 p4(main_cylinder_radius,  -main_cylinder_radius, 0.);
    std::vector<Vec3> bottom_points{p4, p3, p2, p1};
    const auto bottom_surface = simulator.create_point_surface(bottom_points, true);
    fine_cylinder->connect_surface(bottom_surface, 0);
    std::cout << "Normal of bottom surface is " << bottom_surface->get_normal() << "\n";

    // Creating The top plate surface
    Vec3 p5(-main_cylinder_radius, -main_cylinder_radius, filling_height);
    Vec3 p6(-main_cylinder_radius,  main_cylinder_radius, filling_height);
    Vec3 p7(main_cylinder_radius,   main_cylinder_radius, filling_height);
    Vec3 p8(main_cylinder_radius,  -main_cylinder_radius, filling_height);
    std::vector<Vec3> top_points{p5, p6, p7, p8};
    auto top_surface = simulator.create_point_surface(top_points, true);
    fine_cylinder->connect_surface(top_surface, 24);
    std::cout << "Normal of top surface is " << top_surface->get_normal() << "\n";
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
    EngineType::RunForTime run_for_time01(simulator, 1s);
    simulator.run(run_for_time01);
    EngineType::ParticleVelocityLess max_velocity(simulator, 0.1, 0.01s);
    simulator.run(max_velocity);
    simulator.remove_output(output_filling);


    double z_max = 0.055;
    top_surface->move(Vec3(0, 0, z_max - filling_height), Vec3(0, 0, 0));

    top_surface->set_velocity(Vec3(0 ,0, -compaction_velocity));

    auto output_compaction = simulator.create_output(output_directory + "/compaction/", 0.001s);
    output_compaction->print_kinetic_energy = true;
    output_compaction->print_surface_forces = true;
    output_compaction->print_surface_positions = true;
    output_compaction->print_fine_cylinder = true;


    double comp_force = 7.5e6*3.1415*main_cylinder_radius*main_cylinder_radius;
    EngineType::SurfaceNormalForceGreater surface_normal_force_greater(top_surface, comp_force);
    simulator.run(surface_normal_force_greater);
}