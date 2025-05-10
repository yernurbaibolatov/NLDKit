import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science', 'grid'])

script_dir = os.path.dirname(os.path.abspath(__file__))
csv_file = os.path.join(script_dir, '../data/damped_oscillator.csv')


# Load the data
df = pd.read_csv(csv_file)

# Extract parameters for analytic solution
# Assuming the CSV header: time,state_0 (position),state_1 (velocity)
times = df['time'].values
# Define constants (adjust if needed)
omega_0 = 1.0  # natural frequence
gamma = 0.1  # damping constant

# Initial conditions
x0 = df['state_0'].iloc[0]
v0 = df['state_1'].iloc[0]

# Calculate damped natural frequency
omega_d = np.sqrt(max(omega_0**2 - gamma**2, 0))
psi = np.atan2(omega_d*x0, v0+gamma*x0)
A = np.sqrt(x0**2 + ((v0+gamma*x0)/omega_d)**2)
phi = np.atan2(gamma, omega_d)

# Analytic solution: underdamped case
x_analytic = A*np.exp(-gamma * times) * np.sin(omega_d * times + psi)
v_analytic = omega_0*A*np.exp(-gamma*times) * np.cos(omega_d*times + psi + phi)

fig, axes = plt.subplots(2, 2, figsize=(12, 8))

# Position vs time
axes[0, 0].plot(df['time'], df['state_0'], label='Position (Numerical)')
axes[0, 0].plot(df['time'], x_analytic, '--', label='Position (Analytic)')
axes[0, 0].set_xlabel('Time')
axes[0, 0].set_ylabel('Position')
axes[0, 0].set_title('Position vs Time')
axes[0, 0].legend()
axes[0, 0].grid(True)

# Velocity vs time
axes[0, 1].plot(df['time'], df['state_1'], label='Velocity (Numerical)')
axes[0, 1].plot(df['time'], v_analytic, '--', label='Velocity (Analytic)')
axes[0, 1].set_xlabel('Time')
axes[0, 1].set_ylabel('Velocity')
axes[0, 1].set_title('Velocity vs Time')
axes[0, 1].legend()
axes[0, 1].grid(True)

# Phase portrait
axes[1, 0].plot(df['state_0'], df['state_1'])
axes[1, 0].plot(x_analytic, v_analytic, '--')
axes[1, 0].set_xlabel('Position')
axes[1, 0].set_ylabel('Velocity')
axes[1, 0].set_title('Phase Portrait')
axes[1, 0].grid(True)

# Hide the unused subplot
fig.delaxes(axes[1, 1])

plt.tight_layout()
plt.show()