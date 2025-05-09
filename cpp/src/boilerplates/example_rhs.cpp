#include "UserDefinedSystem.hpp"
#include <iostream>

class ExampleSystem : public UserDefinedSystem {
public:
    void rhs(double t, const Eigen::VectorXd& y, Eigen::VectorXd& dydt) override {
        // Simple harmonic oscillator: y'' + y = 0
        dydt[0] = y[1];
        dydt[1] = -y[0];
    }
};

// Expose an instance
ExampleSystem concreteSystem;
DynamicalSystem& mySystem = concreteSystem;