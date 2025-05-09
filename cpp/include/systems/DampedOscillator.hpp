#pragma once

#include "DynamicalSystem.hpp"

class DampedOscillator : public DynamicalSystem {
    double omega_;
    double gamma_;

public:
    // Constructor initializes parameters
    DampedOscillator(double omega, double gamma) {
        setParameter("omega", omega);
        setParameter("gamma", gamma);
    }

    // Right-hand side of the damped oscillator ODE
    void rhs(double t, const Vec& y, Vec& dydt) override {
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
