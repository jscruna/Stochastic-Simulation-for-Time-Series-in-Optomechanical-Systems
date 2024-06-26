import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
csv_data1 = pd.read_csv('variances.csv')

# Extract the data (assuming 'variances' is the column name in the CSV file)
variances = csv_data1.values  # Adjust 'variances5' if the column name is different

# Define constants
T = 10.1
fs = 1e6 / 100
dt = 1 / fs
time = np.arange(0, T, dt)

mass_mean = 7.71e-6  # mean value
kappa_p_mean = 2 * np.pi * 1.64 * 1e6  # mean value
g_mean = -2 * np.pi * 3.2 * 1e4  # mean value
omega_m_mean = 2 * np.pi * 280  # mean value
Detune_mean = 0.0292  # mean value
Temp_mean = 11e-3
n_th_mean = 8e5
gamma_m = 1.1 * 2 * np.pi
N_th = 19

n_x_mean = 2 * gamma_m * (2 * n_th_mean + 1) + 16 * g_mean**2 * (2 * N_th + 1) / ((1 + 4 * Detune_mean**2) * kappa_p_mean)

# Define the decay function
def decay(time):
    return n_x_mean / (2 * gamma_m) * (1 - np.exp(-gamma_m * time))

# Plotting
plt.figure(1)
plt.plot(time, variances, '-', label='Simulation')
plt.plot(time, decay(time), '--', label='Theory')
plt.xlim([0, 10])
plt.xlabel('time (s)')
plt.ylabel('$V_{qq}$ (a.u.)')
plt.legend()
plt.grid(True)
plt.show()


# Load the CSV file
csv_data2 = pd.read_csv('time_series_sample.csv')
trajectories = csv_data2.values

# Plotting
plt.figure(2)
plt.plot(time[1:1000], trajectories[1:1000,0], '-', label='Simulation 1')
plt.plot(time[1:1000], trajectories[1:1000,1], '-', label='Simulation 2')
plt.plot(time[1:1000], trajectories[1:1000,2], '-', label='Simulation 3')
#plt.xlim([0, 10])
plt.xlabel('time (s)')
plt.ylabel('Position $q(t)$ (a.u.)')
plt.legend()
plt.grid(True)
plt.show()
