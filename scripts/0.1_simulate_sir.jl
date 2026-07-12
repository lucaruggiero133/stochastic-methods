using Pkg
Pkg.activate(".")

# Explicitly include the local module file
include("../src/StochasticMethods.jl")
using .StochasticMethods # The dot is important!

using Plots

# --- Parameters ---
N = 1000
avg_degree = 10
λ_real = 0.4
μ_real = 0.2
T = 10

# --- 1. Graph Generation ---
println("Generating sparse graph...")
G = generate_sparse_graph(N, avg_degree)

# --- 2. Simulation ---
println("Simulating ground truth epidemic...")
vT = 1:(T/20):T
x = fill(0, N)
stats = DeltaRecorder(vT)

# Run the Gillespie algorithm
obs_I, obs_R = sample!(x, G, λ_real, μ_real, T; stats=stats)

# Reconstruct absolute numbers via cumulative sum
Ostar_I = cumsum(obs_I)
Ostar_R = cumsum(obs_R)

# --- 3. Plotting ---
println("Plotting results...")
p = scatter(vT, Ostar_I, label="Infected (I)", color=:red, xlabel="Time", ylabel="Individuals", title="SIR Stochastic Trajectory")
scatter!(p, vT, Ostar_R, label="Recovered (R)", color=:blue)

display(p)
println("Done! Close the plot window to exit.")
readline() # Keeps the plot window open until you press Enter