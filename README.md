# Stochastic SIR Epidemics on Sparse Networks

Course project for *Stochastic Simulation Methods* (Politecnico di Torino), written in Julia.

The project simulates an SIR (susceptible–infected–recovered) epidemic as a continuous-time Markov process on a sparse random graph, instead of the mean-field differential equations, and then tries the inverse problem: recover the infection and recovery rates from a noisy observed trajectory.

## Model

Each node of the graph is an individual. Transitions are independent Poisson processes:

- infection, at rate $\lambda$ for each infected neighbour: $S + I \xrightarrow{\lambda} 2I$
- recovery, at rate $\mu$: $I \xrightarrow{\mu} R$

Trajectories are generated exactly with the Gillespie algorithm, using `ExponentialQueues.jl` to sample the next event.

## Inference

Given an observed trajectory $O^\ast$, the likelihood of $(\lambda, \mu)$ is estimated from $M$ simulated trajectories $X_m$. An exact match between discrete trajectories is almost never observed, so the indicator of a match is replaced by a soft weight:

$$\mathbb{P}(O^\ast \mid \lambda, \mu) \approx \frac{1}{M} \sum_{m=1}^{M} \exp\left( -\beta \, \| O^\ast - O(X_m) \| \right),$$

where $\| \cdot \|$ is the Euclidean distance and $\beta$ sets how strict the comparison is. The code works with the log-likelihood, since the likelihood itself underflows over long horizons, and maximises it on a grid of $(\lambda, \mu)$.

## Repository

- `src/` — the module `StochasticMethods` (`Epidemic.jl` for the simulation, `Inference.jl` for the likelihood)
- `scripts/` — the two entry points below
- `notebooks/` — exploration and figures

## Running it

```bash
git clone https://github.com/lucaruggiero133/stochastic-methods.git
cd stochastic-methods
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/0.1_simulate_sir.jl    # simulate and plot one epidemic
julia --project=. scripts/0.2_run_inference.jl   # grid search for (λ, μ)
```
