#include "Integrator.hpp"
#include "Definitions.hpp"
#include <iostream>
#include <cmath>

namespace {
void runIntegration(AbstractDynamicalSystem& system, double dt, Vec& y,
                    double t0, double tf, double t_transient, double output_interval,
                    Mat& results, Vec& times) {
    Vec k1(y.size()), k2(y.size()), k3(y.size()), k4(y.size());
    double t = t0;
    int dim = y.size();
    int num_samples = (output_interval > 0.0)
        ? static_cast<int>(std::floor((tf - t_transient) / output_interval)) + 1
        : 0;
    if (num_samples > 0) {
        results.resize(dim, num_samples);
        times.resize(num_samples);
    }
    int sample_idx = 0;
    double next_output_time = t_transient;

    // Transient phase
    while (t < t_transient) {
        system.rhs(t, y, k1);
        system.rhs(t + 0.5 * dt, y + 0.5 * dt * k1, k2);
        system.rhs(t + 0.5 * dt, y + 0.5 * dt * k2, k3);
        system.rhs(t + dt, y + dt * k3, k4);
        y += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        t += dt;
    }

    // Main phase
    while (t < tf) {
        if (output_interval > 0.0 && t >= next_output_time && sample_idx < num_samples) {
            results.col(sample_idx) = y;
            times(sample_idx) = t;
            next_output_time += output_interval;
            sample_idx++;
        }
        system.rhs(t, y, k1);
        system.rhs(t + 0.5 * dt, y + 0.5 * dt * k1, k2);
        system.rhs(t + 0.5 * dt, y + 0.5 * dt * k2, k3);
        system.rhs(t + dt, y + dt * k3, k4);
        y += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        t += dt;
    }
    std::cout << "Integration complete. Final state:\n" << y.transpose() << std::endl;
}
}

IntegratorDynamic::IntegratorDynamic(AbstractDynamicalSystem& system, double dt)
    : system_(system), dt_(dt) {}

void IntegratorDynamic::setTransientTime(double t_transient) {
    t_transient_ = t_transient;
}

void IntegratorDynamic::setOutputInterval(double interval) {
    output_interval_ = interval;
}

const Mat& IntegratorDynamic::getResults() const {
    return results_;
}

const Vec& IntegratorDynamic::getTimes() const {
    return times_;
}

void IntegratorDynamic::integrate(Vec& y, double t0, double tf) {
    runIntegration(system_, dt_, y, t0, tf, t_transient_, output_interval_, results_, times_);
}

// Integrator2D implementations

namespace {
void runIntegration2D(AbstractDynamicalSystem2D& system, double dt, Vec2& y,
                      double t0, double tf, double t_transient, double output_interval,
                      Mat& results, Vec& times) {
    Vec2 k1, k2, k3, k4;
    double t = t0;
    int num_samples = (output_interval > 0.0)
        ? static_cast<int>(std::floor((tf - t_transient) / output_interval)) + 1
        : 0;
    if (num_samples > 0) {
        results.resize(2, num_samples);
        times.resize(num_samples);
    }
    int sample_idx = 0;
    double next_output_time = t_transient;

    // Transient phase
    while (t < t_transient) {
        system.rhs(t, y, k1);
        system.rhs(t + 0.5 * dt, y + 0.5 * dt * k1, k2);
        system.rhs(t + 0.5 * dt, y + 0.5 * dt * k2, k3);
        system.rhs(t + dt, y + dt * k3, k4);
        y += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        t += dt;
    }

    // Main phase
    while (t < tf) {
        if (output_interval > 0.0 && t >= next_output_time && sample_idx < num_samples) {
            results.col(sample_idx) = y;
            times(sample_idx) = t;
            next_output_time += output_interval;
            sample_idx++;
        }
        system.rhs(t, y, k1);
        system.rhs(t + 0.5 * dt, y + 0.5 * dt * k1, k2);
        system.rhs(t + 0.5 * dt, y + 0.5 * dt * k2, k3);
        system.rhs(t + dt, y + dt * k3, k4);
        y += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        t += dt;
    }
    std::cout << "2D Integration complete. Final state:\n" << y.transpose() << std::endl;
}
}

Integrator2D::Integrator2D(AbstractDynamicalSystem2D& system, double dt)
    : system_(system), dt_(dt) {}

void Integrator2D::setTransientTime(double t_transient) {
    t_transient_ = t_transient;
}

void Integrator2D::setOutputInterval(double interval) {
    output_interval_ = interval;
}

const Mat& Integrator2D::getResults() const {
    return results_;
}

const Vec& Integrator2D::getTimes() const {
    return times_;
}

void Integrator2D::integrate(Vec2& y, double t0, double tf) {
    runIntegration2D(system_, dt_, y, t0, tf, t_transient_, output_interval_, results_, times_);
}

// Integrator3D implementations

namespace {
void runIntegration3D(AbstractDynamicalSystem3D& system, double dt, Vec3& y,
                      double t0, double tf, double t_transient, double output_interval,
                      Mat& results, Vec& times) {
    Vec3 k1, k2, k3, k4;
    double t = t0;
    int num_samples = (output_interval > 0.0)
        ? static_cast<int>(std::floor((tf - t_transient) / output_interval)) + 1
        : 0;
    if (num_samples > 0) {
        results.resize(3, num_samples);
        times.resize(num_samples);
    }
    int sample_idx = 0;
    double next_output_time = t_transient;

    // Transient phase
    while (t < t_transient) {
        system.rhs(t, y, k1);
        system.rhs(t + 0.5 * dt, y + 0.5 * dt * k1, k2);
        system.rhs(t + 0.5 * dt, y + 0.5 * dt * k2, k3);
        system.rhs(t + dt, y + dt * k3, k4);
        y += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        t += dt;
    }

    // Main phase
    while (t < tf) {
        if (output_interval > 0.0 && t >= next_output_time && sample_idx < num_samples) {
            results.col(sample_idx) = y;
            times(sample_idx) = t;
            next_output_time += output_interval;
            sample_idx++;
        }
        system.rhs(t, y, k1);
        system.rhs(t + 0.5 * dt, y + 0.5 * dt * k1, k2);
        system.rhs(t + 0.5 * dt, y + 0.5 * dt * k2, k3);
        system.rhs(t + dt, y + dt * k3, k4);
        y += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        t += dt;
    }
    std::cout << "3D Integration complete. Final state:\n" << y.transpose() << std::endl;
}
}

Integrator3D::Integrator3D(AbstractDynamicalSystem3D& system, double dt)
    : system_(system), dt_(dt) {}

void Integrator3D::setTransientTime(double t_transient) {
    t_transient_ = t_transient;
}

void Integrator3D::setOutputInterval(double interval) {
    output_interval_ = interval;
}

const Mat& Integrator3D::getResults() const {
    return results_;
}

const Vec& Integrator3D::getTimes() const {
    return times_;
}

void Integrator3D::integrate(Vec3& y, double t0, double tf) {
    runIntegration3D(system_, dt_, y, t0, tf, t_transient_, output_interval_, results_, times_);
}

// Integrator wrapper implementations

Integrator::Integrator(AbstractDynamicalSystem& system, double dt)
    : integratorDynamic_(new IntegratorDynamic(system, dt)),
      integrator2D_(nullptr),
      integrator3D_(nullptr) {}

Integrator::Integrator(AbstractDynamicalSystem2D& system, double dt)
    : integratorDynamic_(nullptr),
      integrator2D_(new Integrator2D(system, dt)),
      integrator3D_(nullptr) {}

Integrator::Integrator(AbstractDynamicalSystem3D& system, double dt)
    : integratorDynamic_(nullptr),
      integrator2D_(nullptr),
      integrator3D_(new Integrator3D(system, dt)) {}

void Integrator::setTransientTime(double t_transient) {
    if (integratorDynamic_) {
        integratorDynamic_->setTransientTime(t_transient);
    } else if (integrator2D_) {
        integrator2D_->setTransientTime(t_transient);
    } else if (integrator3D_) {
        integrator3D_->setTransientTime(t_transient);
    }
}

void Integrator::setOutputInterval(double interval) {
    if (integratorDynamic_) {
        integratorDynamic_->setOutputInterval(interval);
    } else if (integrator2D_) {
        integrator2D_->setOutputInterval(interval);
    } else if (integrator3D_) {
        integrator3D_->setOutputInterval(interval);
    }
}

const Mat& Integrator::getResults() const {
    if (integratorDynamic_) {
        return integratorDynamic_->getResults();
    } else {
        static Mat empty_mat;
        return empty_mat;
    }
}


const Vec& Integrator::getTimes() const {
    if (integratorDynamic_) {
        return integratorDynamic_->getTimes();
    } else if (integrator2D_) {
        return integrator2D_->getTimes();
    } else if (integrator3D_) {
        return integrator3D_->getTimes();
    } else {
        static Vec empty_vec;
        return empty_vec;
    }
}

void Integrator::integrate(Vec& y, double t0, double tf) {
    if (integratorDynamic_) {
        integratorDynamic_->integrate(y, t0, tf);
    }
}

void Integrator::integrate(Vec2& y, double t0, double tf) {
    if (integrator2D_) {
        integrator2D_->integrate(y, t0, tf);
    }
}

void Integrator::integrate(Vec3& y, double t0, double tf) {
    if (integrator3D_) {
        integrator3D_->integrate(y, t0, tf);
    }
}