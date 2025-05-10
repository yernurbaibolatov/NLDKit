#include "Integrator.hpp"
#include "Definitions.hpp"
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <vector>

Integrator::Integrator(AbstractDynamicalSystem& system, double dt)
    : system_(system), dt_(dt) {}

void Integrator::setTransientTime(double t_transient) {
    t_transient_ = t_transient;
}

void Integrator::setOutputInterval(double interval) {
    output_interval_ = interval;
}

const Mat& Integrator::getResults() const {
    return results_;
}

const Vec& Integrator::getTimes() const {
    return times_;
}

void Integrator::integrate(Vec& y, double t0, double tf) {
    if (system_.dim == 0) {
        throw std::invalid_argument("Integrator::integrate: system dimension (dim) must be set and positive.");
    }
    if (y.size() != system_.dim) {
        throw std::invalid_argument("Integrator::integrate: input vector y size does not match system dimension.");
    }
    if (k1_.size() != system_.dim) {
        k1_.resize(system_.dim);
        k2_.resize(system_.dim);
        k3_.resize(system_.dim);
        k4_.resize(system_.dim);
    }
    if (y_temp_.size() != system_.dim) {
        y_temp_.resize(system_.dim);
    }
    double t = t0;
    const double dt_half = dt_ / 2.0;
    const double dt_sixth = dt_ / 6.0;
    int dim = system_.dim;
    int num_samples = (output_interval_ > 0.0)
        ? static_cast<int>(std::floor((tf - t_transient_) / output_interval_)) + 1
        : 0;

    int steps_per_sample = 0;
    if (output_interval_ > 0.0) {
        steps_per_sample = static_cast<int>(std::round(output_interval_ / dt_));
        if (steps_per_sample <= 0) steps_per_sample = 1;
    }

    if (num_samples > 0) {
        results_.resize(dim, num_samples);
        times_.resize(num_samples);
    }
    int sample_idx = 0;

    // Transient phase
    while (t < t_transient_) {
        system_.rhs(t, y, k1_);
        y_temp_.noalias() = y;
        y_temp_.noalias() += dt_half * k1_;
        system_.rhs(t + dt_half, y_temp_, k2_);
        y_temp_.noalias() = y;
        y_temp_.noalias() += dt_half * k2_;
        system_.rhs(t + dt_half, y_temp_, k3_);
        y_temp_.noalias() = y;
        y_temp_.noalias() += dt_ * k3_;
        system_.rhs(t + dt_, y_temp_, k4_);
        y.noalias() += dt_sixth * (k1_ + 2.0 * k2_ + 2.0 * k3_ + k4_);
        t += dt_;
    }

    int max_steps = static_cast<int>(std::ceil((tf - t) / dt_));

    // Main phase
    while (sample_idx < num_samples) {
        for (int step = 0; step < steps_per_sample && max_steps > 0; ++step, --max_steps) {
            system_.rhs(t, y, k1_);
            y_temp_.noalias() = y;
            y_temp_.noalias() += dt_half * k1_;
            system_.rhs(t + dt_half, y_temp_, k2_);
            y_temp_.noalias() = y;
            y_temp_.noalias() += dt_half * k2_;
            system_.rhs(t + dt_half, y_temp_, k3_);
            y_temp_.noalias() = y;
            y_temp_.noalias() += dt_ * k3_;
            system_.rhs(t + dt_, y_temp_, k4_);
            y.noalias() += dt_sixth * (k1_ + 2.0 * k2_ + 2.0 * k3_ + k4_);
            t += dt_;
        }
        results_.col(sample_idx) = y;
        times_(sample_idx) = t;
        sample_idx++;
        if (max_steps <= 0) break;
    }
    // std::cout << "Integration complete. Final state:\n" << y.transpose() << std::endl;
}