# Lorenz System Simulations

This repository explores the dynamics of the **Lorenz system**, including classical chaotic behavior and more complex non-stationary scenarios where system parameters change over time.

The goal is to generate reliable datasets and visualizations for studying **chaotic transitions, parameter changes**, and **state forecasting** using both traditional and modern data-driven approaches.

---

## 📘 About the Project

The Lorenz system is a set of three nonlinear differential equations originally introduced to model atmospheric convection.  
Small changes in initial conditions or parameters lead to large differences in trajectories — a hallmark of chaotic behavior.

This project builds on that foundation by adding:
- **Randomized initial conditions** for diverse trajectory sets.  
- **Non-stationary and time-dependent dynamics**, where the Rayleigh number `ρ` changes continuously or stepwise.  
- **Noise injection** to simulate real-world uncertainty.  
- **Flexible simulation modes** for both analytical and machine-learning-based studies.

All simulations output structured data that can be saved to CSV for later analysis or modeling.

---

## ⚙️ Main Components

### 1. `LorenzSystemSim` Class

A configurable simulator for generating trajectories of the Lorenz system under various experimental setups.

**Key features:**
- Adjustable parameters: `σ` (Prandtl number), `β`, and `ρ` (Rayleigh number).  
- Time-dependent `ρ` options (exponential, sinusoidal, or stepwise).  
- Randomized initial conditions and noise control.  
- Multiple run modes for large-scale data generation.  
- CSV export of results for further analysis.

Each simulation produces a structured dataset with columns:
```
SimulationID, SubSimulationID, Time, X, Y, Z, (optional) Rho
```

---

### 2. Simulation Runner Script

Example workflow for generating a large collection of simulations:
- Defines system parameters and initial conditions.  
- Selects simulation types (e.g., increasing/decreasing `ρ`, random order, variable segment lengths).  
- Runs multiple simulations in parallel.  
- Stores all results in a single `.mat` file for convenience.

This script demonstrates how to explore both **stationary** and **non-stationary** regimes of the Lorenz attractor.

---

## 📊 Example Scenarios

| Mode | Description | Use Case |
|------|--------------|----------|
| Simple Run | Classical Lorenz attractor | Visualization or baseline testing |
| Random IC | Randomized initial conditions | Dataset diversity |
| Time-Dependent ρ | Slowly changing system dynamics | Non-stationary behavior |
| Step Change ρ | Abrupt regime transitions | Fault or regime detection |
| Continuous ρ Intervals | Smooth changes across ranges | Long-term forecasting studies |

---

## 📁 Output Format

All simulation outputs can be written to CSV and contain:
- Time series of `x`, `y`, `z` values.
- Corresponding simulation identifiers.
- Optional columns for the varying parameter `ρ`.

These can be used directly in machine-learning pipelines or for numerical analysis (e.g., Lyapunov exponent estimation, bifurcation plots, anomaly detection).

---

## 🔗 References

- [Lorenz System – Wikipedia](https://en.wikipedia.org/wiki/Lorenz_system)  
- [Lorenz (1963) – Deterministic Nonperiodic Flow](https://journals.ametsoc.org/view/journals/atsc/20/2/1520-0469_1963_020_0130_dnf_2_0_co_2.xml)  
- [Tau Simulations Paper (2022)](https://arxiv.org/pdf/2207.00521.pdf)  
- [Using Machine Learning to Predict Statistical Properties of Non-stationary Dynamical Processes (AIP, 2021)](https://pubs.aip.org/aip/cha/article/31/3/033149/342213/Using-machine-learning-to-predict-statistical)

---

## ✨ Notes

- The simulation data can be used for **chaotic time-series forecasting**, **regime detection**, or **state change identification**.  
- For extended studies, the data can be analyzed in Python, MATLAB, or any other environment supporting CSV/Mat files.  
- The class and example script can be adapted to other nonlinear systems by modifying the differential equations.

---

**Author:** Roland Bolboaca  
**Purpose:** Research and educational use — exploring nonlinear and non-stationary dynamical systems.
