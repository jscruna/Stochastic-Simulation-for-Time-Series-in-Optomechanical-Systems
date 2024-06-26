# Stochastic-Simulation-for-Time-Series-in-Optomechanical-Systems
In this project, I developed and implemented a stochastic simulation model to analyze time series data within optomechanical systems. Optomechanical systems, which integrate optical and mechanical components to manipulate and control light and mechanical motion, often exhibit complex behaviors influenced by random fluctuations and noise.
### How to execute it: From your cmd, write the following: C:\Your\File\Address>julia simulation_time_series.jl & python plot.py

In the project, the simulation of the information obtained for the photodetector of a quantum optomechanical system is considered. We validate the information of the optomechanical system through the theoretical relationship obtained from the theoretical value of many samples to calculate the final variance value. In this way, we obtain the variance value of the position as follows:

![Variance Formula](https://latex.codecogs.com/png.latex?V_{qq}(t)=V_{pp}(t)=\frac{n_x}{2\gamma_m}(1-e^{-\gamma_mt}))

Where ![Variance Formula 1](https://latex.codecogs.com/png.latex?V_{qq}(t)) and ![Variance Formula 2](https://latex.codecogs.com/png.latex?V_{pp}(t)) represent the unconditional value of the position and momentum. Additionally, $n_x$ is the Brownian force due to its coupling with light and the thermal bath. $\gamma_m$ is the mechanical dissipation and $t$ is time. The solution to this equation shows how the variance of a set of samples changes over time.

## Plotting Results

The plotting result is in the file "Results" from this repository

The numerical simulation is based on the paper: Santiago-Condori, J. G., Yamamoto, N., Matsumoto, N.(2023). Verification of conditional mechanical squeezing for a mg-scale pendulum near quantum regimes (arXiv:2008.10848).
