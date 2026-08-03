//
// Created by erolsson on 8/3/26.
//

#include "stone_material.h"

#include "../utilities/file_reading_functions.h"
#include "../utilities/printing_functions.h"

using namespace DEM;
StoneMaterial::StoneMaterial(const ParameterMap& parameters) :
    MaterialBase(parameters),
    E(parameters.get_parameter("E")),
    nu(parameters.get_parameter("nu")),
    unloading_exponent(parameters.get_parameter("unloading_exponent")),
    mu(parameters.get_parameter("mu")),
    mu_wall(parameters.get_parameter("mu_wall")),
    hs(parameters.get_parameter("hs")),
    Fs(parameters.get_parameter("Fs")),
    min_crack_distance(parameters.get_parameter("min_crack_distance")),
    weibull_fracture_stress(parameters.get_parameter("weibull_fracture_stress")),
    weibull_exponent(parameters.get_parameter("weibull_exponent")),
    weibull_ref_volume(parameters.get_parameter("weibull_ref_volume")),
    plastic_linear_depth(parameters.get_parameter("plastic_linear_depth"))
{

}

std::string StoneMaterial::restart_data() const {
    std::ostringstream ss;
    ss << named_print("stone_material", "type") << ", "
       << MaterialBase::restart_data() << ", "
       << named_print(E, "E") << ", "
       << named_print(nu, "nu") << ", "
       << named_print(unloading_exponent, "unloading_exponent") << ", "
       << named_print(mu, "mu") << ", "
       << named_print(mu_wall, "mu_wall") << ", "
       << named_print(hs, "hs") << ", "
       << named_print(Fs, "Fs") << ", "
       << named_print(min_crack_distance, "min_crack_distance") << ", "
       << named_print(weibull_fracture_stress, "weibull_fracture_stress") << ", "
       << named_print(weibull_exponent, "weibull_exponent") << ", "
       << named_print(weibull_ref_volume, "weibull_ref_volume") << ", "
       << named_print(plastic_linear_depth, "plastic_linear_depth");
    return ss.str();
}