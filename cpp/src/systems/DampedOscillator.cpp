#include "DampedOscillator.hpp"

// Create an instance of the damped oscillator with desired parameters
DampedOscillator damped(omega = 1.0, gamma = 0.1);

// Expose a reference to the system for integration
AbstractDynamicalSystem& mySystem = damped;