#pragma once

#include "AbstractDynamicalSystem.hpp"
#include <map>
#include <string>
#include <stdexcept>

class DynamicalSystem : public AbstractDynamicalSystem {
public:
    DynamicalSystem() = default;
    ~DynamicalSystem() override = default;

    void setParameter(const std::string& name, double value) override {
        parameters_[name] = value;
    }

    double getParameter(const std::string& name) const override {
        auto it = parameters_.find(name);
        if (it != parameters_.end()) {
            return it->second;
        } else {
            throw std::invalid_argument("Unknown parameter: " + name);
        }
    }

    /**
     * @brief Computes the right-hand side of the dynamical system.
     * 
     * IMPORTANT: The sizes of the input vectors y and dydt must exactly match the system's dimension (dim).
     * It is the responsibility of the concrete system implementation to ensure that y.size() == dim 
     * and dydt.size() == dim. Any mismatch may result in undefined behavior or runtime errors.
     *
     * @param t Current time
     * @param y State vector (input), must have size equal to dim
     * @param dydt Derivative vector (output), must have size equal to dim
     */
    void rhs(double t, const Vec& y, Vec& dydt) override = 0;

protected:
    std::map<std::string, double> parameters_;
};