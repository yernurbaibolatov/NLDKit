#include <iostream>
#include "Definitions.hpp"
#include "Integrator.hpp"
#include "systems/DampedOscillator.hpp"

int main() {
    // Define system parameters
    double omega = 1.0;   // stiffness (natural frequency)
    double gamma = 0.1;   // damping coefficient

    // Create the damped oscillator system
    DampedOscillator oscillator(omega, gamma);

    // Set initial conditions: position = 1.0, velocity = 0.0
    Vec y(2);
    y(0) = 0.2;  // x
    y(1) = -2.0;  // v

    // Set up the integrator
    double dt = 0.001;
    double t0 = 0.0;
    double tf = 50.0;

    Integrator integrator(oscillator, dt);
    integrator.setTransientTime(0.0);
    integrator.setOutputInterval(0.1);

    // Run integration
    integrator.integrate(y, t0, tf);

    // Print the results
    const Mat& results = integrator.getResults();
    const Vec& times = integrator.getTimes();

    std::cout << "Finished integrating" << std::endl;
    integrator.writeResultsToCSV("../../data/damped_oscillator.csv");

    return 0;
}