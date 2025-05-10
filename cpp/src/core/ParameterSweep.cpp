#include "ParameterSweep.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

ParameterSweep::ParameterSweep(AbstractDynamicalSystem& system, const std::string& param_name)
    : system_(system), param_name_(param_name) {}

void ParameterSweep::setParameterRange(double start, double end, int num_steps) {
    param_values_.clear();
    if (num_steps < 2) {
        param_values_.push_back(start);
    } else {
        double step = (end - start) / (num_steps - 1);
        for (int i = 0; i < num_steps; ++i) {
            param_values_.push_back(start + i * step);
        }
    }
}

void ParameterSweep::setTransientTime(double transient_time) {
    transient_time_ = transient_time;
}

void ParameterSweep::setPostProcessingFunction(PostProcessFunc func) {
    post_process_ = func;
}

void ParameterSweep::runSweep(const Vec& y0, double t0, double tf) {
    if (!post_process_) {
        std::cerr << "Post-processing function is not set!" << std::endl;
        return;
    }

    int num_params = static_cast<int>(param_values_.size());
    std::vector<Vec> all_results;

    for (double param_value : param_values_) {
        system_.setParameter(param_name_, param_value);

        Integrator integrator(system_, /* dt = */ 0.01);  // You can expose dt as configurable if needed
        integrator.setTransientTime(transient_time_);

        Vec y_copy = y0;
        integrator.integrate(y_copy, t0, tf);
        const Mat& result = integrator.getResults();

        Vec processed = post_process_(result);
        all_results.push_back(processed);
    }

    // Store as a matrix: each column is the processed result for a parameter value
    if (!all_results.empty()) {
        int result_dim = static_cast<int>(all_results[0].rows());
        processed_results_.resize(result_dim, num_params);
        for (int i = 0; i < num_params; ++i) {
            processed_results_.col(i) = all_results[i];
        }
    }
}

const std::vector<double>& ParameterSweep::getParameterValues() const {
    return param_values_;
}

const Mat& ParameterSweep::getProcessedResults() const {
    return processed_results_;
}

void ParameterSweep::writeResultsToCSV(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    file << std::scientific << std::setprecision(10);

    // Header
    file << param_name_;
    for (int i = 0; i < processed_results_.rows(); ++i) {
        file << ",Result_" << i;
    }
    file << "\n";

    // Data rows
    for (int j = 0; j < processed_results_.cols(); ++j) {
        file << param_values_[j];
        for (int i = 0; i < processed_results_.rows(); ++i) {
            file << "," << processed_results_(i, j);
        }
        file << "\n";
    }

    file.close();
}