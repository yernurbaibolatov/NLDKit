#pragma once

#include <Eigen/Dense>

class DynamicalSystem {
public:
    virtual ~DynamicalSystem() = default;
    virtual void rhs(double t, const Eigen::VectorXd& y, Eigen::VectorXd& dydt) = 0;
};