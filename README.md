This repository implements a discrete stochastic simulation of the Susceptible-Infected-Recovered (SIR) epidemic model on sparse graphs, moving beyond traditional mean-field Ordinary Differential Equations (ODEs) to capture realistic, event-driven contagion dynamics.It features a custom physics engine built in Julia using the Gillespie Algorithm (Stochastic Simulation Algorithm), alongside an inference engine that utilizes Monte Carlo Maximum Log-Likelihood Estimation to reverse-engineer hidden transmission parameters from noisy, discrete observations.📌 FeaturesExact Stochastic Simulation: Continuous-time Markov Process modeling using ExponentialQueues for highly efficient event sampling.Network Dynamics: Epidemic propagation over randomly generated sparse boolean networks (preventing self-infection and limiting contact degrees).Parameter Inference: A custom Monte Carlo Maximum Likelihood grid-search algorithm to estimate the true infection and recovery rates from limited trajectory data.Modular Architecture: Production-ready code structure separating the mathematical backend from execution scripts and analytical notebooks.📐 Mathematical Formulation1. The Stochastic SIR ModelIndividuals are represented as nodes on a sparse adjacency matrix. State transitions occur continuously in time as independent Poisson processes:Infection: A susceptible individual interacting with an infected neighbor becomes infected at rate $\\lambda$.$$S + I \\xrightarrow{\\lambda} 2I$$Recovery: An infected individual recovers at a constant rate $\\mu$.$$I \\xrightarrow{\\mu} R$$2. Log-Likelihood Inference via Soft ApproximationTo deduce the true parameters from an observed trajectory $\\mathcal{O}^\*$, we sample $M$ Monte Carlo trajectories $X\_m$ for different $(\\lambda, \\mu)$ pairs. Because finding exact discrete trajectory matches is computationally intractable for large state spaces, we approximate the Kronecker delta using a continuous, distance-based soft weight:$$\\mathbb{P}(O^\*\\vert{}\\lambda,\\mu) \\approx \\frac{1}{M}\\sum\_{m=1}^M e^{-\\beta \\vert{}\\vert{}O^\*-O(X\_m)\\vert{}\\vert{}}$$where $\\vert{}\\vert{}\\cdot\\vert{}\\vert{}$ is the Euclidean norm and $\\beta$ controls the strictness of the metric. To prevent numerical underflow over long time horizons, the inference engine computes the Log-Likelihood.📂 Repository ArchitecturePlaintext├── src/

│   ├── StochasticMethods.jl    # Main module wrapper and exports

│   ├── Epidemic.jl             # Graph generation and Gillespie SSA core

│   └── Inference.jl            # Likelihood and distance metrics

├── scripts/

│   ├── 0.1\_simulate\_sir.jl     # Generates ground truth and visualizes trajectory

│   └── 0.2\_run\_inference.jl    # Executes the Monte Carlo grid search

├── notebooks/

│   └── Stochastic\_SIR\_Inference.ipynb  # Comprehensive narrative and results

├── Project.toml                # Julia environment and dependencies

└── README.md

🚀 Quick Start1. Clone the Repository \& Instantiate the EnvironmentBashgit clone https://github.com/yourusername/stochastic-sir-inference.git

cd stochastic-sir-inference

julia

Inside the Julia REPL, press ] to enter the package manager and instantiate the environment:Juliapkg> activate .

pkg> instantiate

2\. Run the Physical SimulationSimulate a ground truth epidemic trajectory and visualize the time-series:Bashjulia scripts/0.1\_simulate\_sir.jl

3\. Run the Inference EngineExecute the Monte Carlo log-likelihood grid search to estimate $\\lambda^\*$ and $\\mu^\*$:Bashjulia scripts/0.2\_run\_inference.jl

📊 Results(Note: Once you run the scripts and generate the plots, save the Heatmap and Trajectory images into an assets/ folder and link them here to show off the visual results of your inference engine.)

