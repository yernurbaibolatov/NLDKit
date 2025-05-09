#include "Integrator.hpp"
#include <iostream>

Integrator::Integrator(DynamicalSystem& system, double dt)
    : system_(system), dt_(dt) {}

void Integrator::integrate(Eigen::VectorXd& y, double t0, double tf) {
    double t = t0;
    Eigen::VectorXd dydt(y.size());

    while (t < tf) {
        system_.rhs(t, y, dydt);
        y += dt_ * dydt;  // simple Euler (can be replaced with RK4)
        t += dt_;
    }

    std::cout << "Integration complete. Final state:\n" << y << std::endl;
}