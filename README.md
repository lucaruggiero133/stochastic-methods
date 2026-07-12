# Stochastic SIR Modeling and Parameter Inference on Sparse Networks

This repository implements a discrete stochastic simulation of the **Susceptible-Infected-Recovered (SIR)** epidemic model on sparse graphs, moving beyond traditional mean-field Ordinary Differential Equations (ODEs) to capture realistic, event-driven contagion dynamics. 

It features a custom physics engine built in Julia using the **Gillespie Algorithm (Stochastic Simulation Algorithm)**, alongside an inference engine that utilizes Monte Carlo Maximum Log-Likelihood Estimation to reverse-engineer hidden transmission parameters from noisy, discrete observations.

## 📌 Features
* **Exact Stochastic Simulation:** Continuous-time Markov Process modeling using `ExponentialQueues` for highly efficient event sampling.
* **Network Dynamics:** Epidemic propagation over randomly generated sparse boolean networks (preventing self-infection and limiting contact degrees).
* **Parameter Inference:** A custom Monte Carlo Maximum Likelihood grid-search algorithm to estimate the true infection and recovery rates from limited trajectory data.
* **Modular Architecture:** Production-ready code structure separating the mathematical backend from execution scripts and analytical notebooks.

## 📐 Mathematical Formulation

### 1. The Stochastic SIR Model
Individuals are represented as nodes on a sparse adjacency matrix. State transitions occur continuously in time as independent Poisson processes:

* **Infection:** A susceptible individual interacting with an infected neighbor becomes infected at rate $\lambda$.

$$S + I \xrightarrow{\lambda} 2I$$

* **Recovery:** An infected individual recovers at a constant rate $\mu$.

$$I \xrightarrow{\mu} R$$

### 2. Log-Likelihood Inference via Soft Approximation
To deduce the true parameters from an observed trajectory $O^*$, we sample $M$ Monte Carlo trajectories $X_m$ for different $(\lambda, \mu)$ pairs. Because finding exact discrete trajectory matches is computationally intractable for large state spaces, we approximate the Kronecker delta using a continuous, distance-based soft weight:

$$\mathbb{P}(O^* \mid \lambda, \mu) \approx \frac{1}{M} \sum_{m=1}^{M} \exp\left( -\beta \| O^* - O(X_m) \| \right)$$

where $\| \cdot \|$ is the Euclidean norm and $\beta$ controls the strictness of the metric. To prevent numerical underflow over long time horizons, the inference engine computes the **Log-Likelihood**.

## 🚀 Quick Start

**1. Clone the Repository & Instantiate the Environment**
```bash
git clone [https://github.com/yourusername/stochastic-sir-inference.git](https://github.com/yourusername/stochastic-sir-inference.git)
cd stochastic-sir-inference
julia
```

Inside the Julia REPL, press ] to enter the package manager and instantiate the environment:
```bash
pkg> activate .
pkg> instantiate
```

**2. Run the Physical Simulation**
Simulate a ground truth epidemic trajectory and visualize the time-series:
```bash
julia scripts/0.1_simulate_sir.jl
```

**3. Run the Inference Engine**
Execute the Monte Carlo log-likelihood grid search to estimate $\lambda^{*$ and $\mu^{*$:
```bash
julia scripts/0.2_run_inference.jl
```


