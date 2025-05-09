#pragma once

#include "DynamicalSystem.hpp"

class Integrator {
public:
    Integrator(DynamicalSystem& system, double dt);
    void integrate(Eigen::VectorXd& y, double t0, double tf);

private:
    DynamicalSystem& system_;
    double dt_;
};