//
// Created by erolsson on 3/29/26.
//

#ifndef DEMSIM_FINECYLINDER_H
#define DEMSIM_FINECYLINDER_H

#include <vector>

class Particle;
class Surface;
class Cylinder;
class Rectangle;

class EvaluationPoint{
    double z0_ = 0;
    double z1_ = 0;
};


class FineCylinder {
    std::vector<EvaluationPoint> evaluation_points_;
    const std::vector<Particle*> particles_;
    const Cylinder* cylinder_ = nullptr;
    const Rectangle* top_surface_ = nullptr;
    const Rectangle* bottom_surface = nullptr;

    void update();
};






#endif //DEMSIM_FINECYLINDER_H