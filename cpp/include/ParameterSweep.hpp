#pragma once

#include "AbstractDynamicalSystem.hpp"
#include "Integrator.hpp"
#include "Definitions.hpp"

#include <vector>
#include <string>
#include <functional>

class ParameterSweep {
public:
    using PostProcessFunc = std::function<Vec(const Mat& result)>;

    ParameterSweep(AbstractDynamicalSystem& system, const std::string& param_name);

    void setParameterRange(double start, double end, int num_steps);
    void setTimeStep(double dt);
    void setTransientTime(double transient_time);
    void setPostProcessingFunction(PostProcessFunc func);

    void runSweep(const Vec& y0, double t0, double tf);
    const std::vector<double>& getParameterValues() const;
    const Mat& getProcessedResults() const;

    // Save parameter values and processed results to CSV
    void writeResultsToCSV(const std::string& filename) const;

private:
    AbstractDynamicalSystem& system_;
    std::string param_name_;
    std::vector<double> param_values_;
    double transient_time_ = 0.0;
    double dt_ = Constants::DEFAULT_DT;  // Default time step from Definitions.hpp

    PostProcessFunc post_process_;
    Mat processed_results_;  // Each column: processed result for each parameter
};
