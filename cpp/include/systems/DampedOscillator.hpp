#pragma once

#include "DynamicalSystem.hpp"

// Class representing a damped harmonic oscillator system
class DampedOscillator : public DynamicalSystem {
    double omega_;  // Natural frequency
    double gamma_;  // Damping coefficient

public:
    // Constructor initializes system dimension and parameters
    DampedOscillator(double omega, double gamma) {
        dim = 2;
        setParameter("omega", omega);
        setParameter("gamma", gamma);
    }

    // Compute the right-hand side of the ODE
    void rhs(double t, const Vec& y, Vec& dydt) override {
        dydt.resize(2);
        dydt[0] = y[1];
        dydt[1] = -2.0 * gamma_ * y[1] - omega_ * omega_ * y[0];
    }

    // Set parameter by name and update cached variables
    void setParameter(const std::string& name, double value) override {
        DynamicalSystem::setParameter(name, value);
        if (name == "omega") {
            omega_ = value;
        } else if (name == "gamma") {
            gamma_ = value;
        } else {
            throw std::invalid_argument("Unknown parameter: " + name);
        }
    }
};
