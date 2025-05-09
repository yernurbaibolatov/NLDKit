#pragma once

#include "AbstractDynamicalSystem.hpp"
#include "Definitions.hpp"

// IntegratorDynamic class: performs numerical integration using the RK4 method,
// with support for transient time, controlled output sampling, and result storage.

class IntegratorDynamic {
public:
    IntegratorDynamic(AbstractDynamicalSystem& system, double dt);  // Constructor with system reference and time step
    void integrate(Vec& y, double t0, double tf);            // Runs integration from t0 to tf

    void setTransientTime(double t_transient);               // Sets transient phase duration (no recording)
    void setOutputInterval(double interval);                 // Sets output sampling interval
    const Mat& getResults() const;                           // Returns matrix of saved results (dim × num_samples)
    const Vec& getTimes() const;                             // Returns vector of saved timestamps

private:
    AbstractDynamicalSystem& system_;  // Reference to the system being integrated
    double dt_;                        // Integration time step

    double t_transient_ = 0.0;         // Transient integration time
    double output_interval_ = -1.0;    // Output sampling interval; -1 means no output recording
    Mat results_;                      // Stores results: (dim × num_samples)
    Vec times_;                        // Stores corresponding times: (num_samples)
};

// Integrator2D: performs RK4 integration for 2D systems using Vec2
class Integrator2D {
public:
    Integrator2D(AbstractDynamicalSystem2D& system, double dt);  // Constructor with system reference and time step
    void integrate(Vec2& y, double t0, double tf);               // Runs integration from t0 to tf

    void setTransientTime(double t_transient);                   // Sets transient phase duration (no recording)
    void setOutputInterval(double interval);                     // Sets output sampling interval
    const Mat& getResults() const;                               // Returns matrix of saved results (2 × num_samples)
    const Vec& getTimes() const;                                 // Returns vector of saved timestamps

private:
    AbstractDynamicalSystem2D& system_;  // Reference to the 2D system
    double dt_;                           // Integration time step

    double t_transient_ = 0.0;           // Transient integration time
    double output_interval_ = -1.0;      // Output sampling interval; -1 means no output recording
    Mat results_;                        // (2 × num_samples)
    Vec times_;                          // (num_samples)
};

// Integrator3D: performs RK4 integration for 3D systems using Vec3
class Integrator3D {
public:
    Integrator3D(AbstractDynamicalSystem3D& system, double dt);  // Constructor with system reference and time step
    void integrate(Vec3& y, double t0, double tf);               // Runs integration from t0 to tf

    void setTransientTime(double t_transient);                   // Sets transient phase duration (no recording)
    void setOutputInterval(double interval);                     // Sets output sampling interval
    const Mat& getResults() const;                               // Returns matrix of saved results (3 × num_samples)
    const Vec& getTimes() const;                                 // Returns vector of saved timestamps

private:
    AbstractDynamicalSystem3D& system_;  // Reference to the 3D system
    double dt_;                           // Integration time step

    double t_transient_ = 0.0;           // Transient integration time
    double output_interval_ = -1.0;      // Output sampling interval; -1 means no output recording
    Mat results_;                        // (3 × num_samples)
    Vec times_;                          // (num_samples)
};

// Integrator wrapper: automatically selects the correct integrator based on the system type
class Integrator {
public:
    Integrator(AbstractDynamicalSystem& system, double dt);
    Integrator(AbstractDynamicalSystem2D& system, double dt);
    Integrator(AbstractDynamicalSystem3D& system, double dt);

    // Overloaded integrate methods for different vector types
    void integrate(Vec& y, double t0, double tf);
    void integrate(Vec2& y, double t0, double tf);
    void integrate(Vec3& y, double t0, double tf);

    void setTransientTime(double t_transient);
    void setOutputInterval(double interval);

    const Mat& getResults() const;
    const Vec& getTimes() const;

private:
    std::unique_ptr<IntegratorDynamic> integratorDynamic_;
    std::unique_ptr<Integrator2D> integrator2D_;
    std::unique_ptr<Integrator3D> integrator3D_;
};