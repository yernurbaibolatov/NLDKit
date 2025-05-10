#pragma once

#include "AbstractDynamicalSystem.hpp"
#include "Definitions.hpp"

// Integrator class: performs numerical integration using the RK4 method,
// with support for transient time, controlled output sampling, and result storage.

class Integrator {
public:
    Integrator(AbstractDynamicalSystem& system, double dt);  // Constructor with system reference and time step
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

    // Internal RK4 buffers (pre-allocated for performance)
    mutable Vec k1_, k2_, k3_, k4_, y_temp_;
};