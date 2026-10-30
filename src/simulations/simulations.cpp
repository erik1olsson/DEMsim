//
// Created by erolsson on 2018-09-18.
//

#include <map>

#include "simulations.h"

std::map<std::string, DEM::SimulationFunctionPtr> DEM::valid_simulations() {
    // Add all valid simulation routines here
    return std::map<std::string, SimulationFunctionPtr>{
        //{"gyratory_compaction",            gyratory_compaction},
        //    {"closed_die_compaction",          closed_die_compaction},
        //    {"contact_tester",                 contact_tester},
        //    {"cyclic_triaxial",                cyclic_triaxial},
        //    {"proctor_test",                   proctor_test},
        //    {"stone_compaction",               stone_compaction},
        //    {"electrode_compaction",           electrode_compaction},
        //    {"electrode_mechanical_test",      electrode_mechanical_test},
        //    {"electrode_cylinder_compaction",  electrode_cylinder_compaction},
        //    {"periodic_bc_tester",             periodic_bc_tester},
        //    {"periodic_bc_simulation",         periodic_bc_simulation},
        //    {"filling_periodic_box",           filling_periodic_box},
        //    {"battery_rve_compaction",         battery_rve_compaction},
        //    {"restart_electrode",              restart_electrode},
        //    {"porous_electrode_rve",           porous_electrode_rve},
        //    {"asphalt_shear_box",              asphalt_shear_box},
        //    {"asphalt_shear_box_bonded",       asphalt_shear_box_bonded},
            {"kcv_uniaxial",                   kcv_uniaxial},
            {"kcv_fines",                      kcv_fines},
            {"kcv_particle_resistance",        kcv_particle_resistance},
            {"kcv_particle_response",          kcv_particle_response},
            {"kcv_particle_drop",              kcv_particle_drop},
    };
}

