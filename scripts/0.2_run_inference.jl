using Pkg
Pkg.activate(".")

# Explicitly include the local module file
include("../src/StochasticMethods.jl")
using .StochasticMethods # The dot is important!

using Plots

# --- 1. Setup Ground Truth ---
N = 1000
G = generate_sparse_graph(N, 10)
λ_real, μ_real = 0.4, 0.2
T = 10
vT = 1:(T/20):T

x = fill(0, N)
stats = DeltaRecorder(vT)
obs_I, obs_R = sample!(x, G, λ_real, μ_real, T; stats=stats)
Ostar = (cumsum(obs_I), cumsum(obs_R))

# --- 2. Inference Parameters ---
M = 200     # Monte Carlo samples per parameter pair
β = 1.0     # Metric strictness (β = 1/2σ²)

λ_range = 0.1:0.1:1.1
μ_range = 0.1:0.1:1.1

# --- 3. Run Maximum Log-Likelihood ---
println("Running Log-Likelihood Inference over grid (this may take a minute)...")
Llog = maxloglikelihood(λ_range, μ_range, G, T, M, Ostar, β)

# Find the best parameters
max_idx = argmax(Llog)
λ_star = λ_range[max_idx[1]]
μ_star = μ_range[max_idx[2]]

println("\n--- Inference Complete ---")
println("Ground truth: λ = $λ_real, μ = $μ_real")
println("Estimated:    λ* = $λ_star, μ* = $μ_star")

# --- 4. Plot Heatmap ---
p = heatmap(
    λ_range, μ_range, Llog', # Transposed for correct spatial alignment
    title="Log-Likelihood (λ*=$λ_star, μ*=$μ_star)",
    xlabel="λ (Infection Rate)",
    ylabel="μ (Recovery Rate)",
    seriescolor=:thermal
)

display(p)
println("Press Enter to close.")
readline()