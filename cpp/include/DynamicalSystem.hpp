#pragma once

#include "AbstractDynamicalSystem.hpp"
#include <map>
#include <string>
#include <stdexcept>

class ParameterManager {
protected:
    std::map<std::string, double> parameters_;

public:
    virtual ~ParameterManager() = default;

    void setParameter(const std::string& name, double value) {
        parameters_[name] = value;
    }

    double getParameter(const std::string& name) const {
        auto it = parameters_.find(name);
        if (it != parameters_.end()) {
            return it->second;
        } else {
            throw std::invalid_argument("Unknown parameter: " + name);
        }
    }
};

class DynamicalSystem : public AbstractDynamicalSystem, public ParameterManager {
public:
    DynamicalSystem() = default;
    ~DynamicalSystem() override = default;

    void rhs(double t, const Vec& y, Vec& dydt) override = 0;
};

// Dynamical system class for 2D systems
class DynamicalSystem2D : public AbstractDynamicalSystem2D, public ParameterManager {
public:
    DynamicalSystem2D() = default;
    ~DynamicalSystem2D() override = default;

    void rhs(double t, const Vec2& y, Vec2& dydt) override = 0;
};

// Dynamical system class for 3D systems
class DynamicalSystem3D : public AbstractDynamicalSystem3D, public ParameterManager {
public:
    DynamicalSystem3D() = default;
    ~DynamicalSystem3D() override = default;

    void rhs(double t, const Vec3& y, Vec3& dydt) override = 0;
};