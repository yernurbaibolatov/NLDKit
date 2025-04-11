import numpy as np

def rk4_step(f, t, x, h):
    """Performs a single Runge-Kutta 4th order step."""
    k1 = f(t, x)
    k2 = f(t + h/2, x + h/2 * k1)
    k3 = f(t + h/2, x + h/2 * k2)
    k4 = f(t + h, x + h * k3)
    return x + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

def solve_rk4(f, x0, t0=0.0, tmax=50.0, dt=0.01):
    """
    Solves a system of ODEs using the fixed-step Runge-Kutta 4th order method.
    
    Parameters:
        f : function
            Function of the form f(t, x) returning dx/dt.
        x0 : array_like
            Initial state vector.
        t0 : float
            Initial time.
        tmax : float
            Maximum integration time.
        dt : float
            Time step.

    Returns:
        ts : ndarray
            Array of time points.
        xs : ndarray
            Array of state vectors corresponding to each time point.
    """
    ts = [t0]
    xs = [np.array(x0)]
    t = t0
    x = np.array(x0)

    while t < tmax:
        x = rk4_step(f, t, x, dt)
        t += dt
        ts.append(t)
        xs.append(x.copy())

    return np.array(ts), np.array(xs)
