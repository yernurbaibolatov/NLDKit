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

// Specialization for 2D systems
class AbstractDynamicalSystem2D {
public:
    virtual ~AbstractDynamicalSystem2D() = default;

    static constexpr int dim = 2;  // Fixed dimension for 2D systems
    virtual void rhs(double t, const Vec2& y, Vec2& dydt) = 0;

    virtual void setParameter(const std::string& name, double value) {
        throw std::runtime_error("setParameter() not implemented for this 2D system.");
    }

    virtual double getParameter(const std::string& name) const {
        throw std::runtime_error("getParameter() not implemented for this 2D system.");
    }
};

// Specialization for 3D systems
class AbstractDynamicalSystem3D {
public:
    virtual ~AbstractDynamicalSystem3D() = default;

    static constexpr int dim = 3;  // Fixed dimension for 3D systems
    virtual void rhs(double t, const Vec3& y, Vec3& dydt) = 0;

    virtual void setParameter(const std::string& name, double value) {
        throw std::runtime_error("setParameter() not implemented for this 3D system.");
    }

    virtual double getParameter(const std::string& name) const {
        throw std::runtime_error("getParameter() not implemented for this 3D system.");
    }
};