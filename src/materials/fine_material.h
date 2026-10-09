//
// Created by erolsson on 8/17/26.
//

#ifndef DEMSIM_FINEMATERIAL_H
#define DEMSIM_FINEMATERIAL_H

#include <array>

#include "material_base.h"

namespace DEM {
    class FineMaterial : public MaterialBase {
    public:
        FineMaterial(unsigned id_number, double density) :
        MaterialBase(id_number, density)
        {}
        explicit FineMaterial(const ParameterMap& parameters);
        ~FineMaterial() override = default;
        [[nodiscard]] std::string restart_data() const override;

        double db = 0;
        double D_th = 0;
        double kl = 0;
        double n = 0;

        double a = 0;
        double b = 0;

        void assign_parameters(const std::array<double, 6>& parameters) {
            db = parameters[0];
            D_th = parameters[1];
            kl = parameters[2];
            n = parameters[3];

            a = parameters[4];
            b = parameters[5];
        }

    };



    inline std::array<double, 6> fuller5_25_parameters = {
        6.65180822e-04, 5.15141338e-01, 1.31312298e+05*1e6, 3.18855802e+00, 0.658433977547588, 0.5
    };

}



#endif //DEMSIM_FINEMATERIAL_H
