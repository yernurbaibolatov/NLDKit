#pragma once

#include "Definitions.hpp"
#include <string>
#include <stdexcept>

class AbstractDynamicalSystem {
public:

    virtual ~AbstractDynamicalSystem() = default;

    int dim = 0;  // System dimension; must be set by each concrete system.
    virtual void rhs(double t, const Vec& y, Vec& dydt) = 0;

    // Set a named parameter to a value (default: throws if unsupported)
    virtual void setParameter(const std::string& name, double value) {
        throw std::runtime_error("setParameter() not implemented for this system.");
    }

    // Get a named parameter's value (default: throws if unsupported)
    virtual double getParameter(const std::string& name) const {
        throw std::runtime_error("getParameter() not implemented for this system.");
    }
};
