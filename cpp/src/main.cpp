#include "Integrator.hpp"
#include "UserDefinedSystem.hpp"
#include <Eigen/Dense>
#include <iostream>

// Externally defined system (from boilerplate)
extern DynamicalSystem& mySystem;

int main() {
    Eigen::VectorXd y(2);
    y << 1.0, 0.0;  // Initial condition: position=1, velocity=0

    double t0 = 0.0;
    double tf = 10.0;
    double dt = 0.01;

    Integrator integrator(mySystem, dt);
    integrator.integrate(y, t0, tf);

    return 0;
}