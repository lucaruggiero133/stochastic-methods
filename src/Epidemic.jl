# Helper functions
rand_exponential(r) = -log(rand()) / r
neighbors(Λ, i) = @views rowvals(Λ)[nzrange(Λ, i)]

"""
    generate_sparse_graph(N, avg_degree)

Generates a random sparse boolean network with no self-loops.
"""
function generate_sparse_graph(N::Int, avg_degree::Real)
    G = sprand(Bool, N, N, avg_degree / N)
    G -= Diagonal(G) # Ensures no individual can infect themselves
    return sparse(G)
end

"""
    DeltaRecorder

A callable struct that safely tracks the deltas of infected and recovered 
individuals at discrete times defined by `ran`.
"""
mutable struct DeltaRecorder{T<:AbstractVector}
    vI::Vector{Int}
    vR::Vector{Int}
    ran::T
    k::Int
end

# Constructor
function DeltaRecorder(ran::AbstractVector)
    return DeltaRecorder(fill(0, length(ran)), fill(0, length(ran)), ran, 1)
end

# Makes the struct callable like a function: stats(t, x, i)
function (r::DeltaRecorder)(t, x, i)
    while r.k <= length(r.ran) && r.ran[r.k] < t
        r.k += 1
    end
    if r.k <= length(r.ran)
        r.vI[r.k] += (x[i] == 1) - (x[i] == 2)
        r.vR[r.k] += (x[i] == 2)
    end
end

"""
    sample!(x, G, λ, μ, T; stats)

Simulates the SIR epidemic process on a graph G using the Gillespie algorithm (SSA).
"""
function sample!(x, G, λ, μ, T; stats = (t, x, i) -> nothing)
    N = size(G, 1)
    t = 0.0

    # Sample initial state
    x .= 0 
    x[rand(1:N)] = 1 
    
    Q = ExponentialQueue(N) 
    for i in 1:N 
        if x[i] == 1
            x[i] = 0
            Q[i] = Inf 
        end
    end

    while !isempty(Q)
        i, Δt = pop!(Q) 
        t += Δt
        t > T && break
        
        x[i] += 1
        if x[i] == 1 # Infection event
            for j in neighbors(G, i)
                if x[j] == 0
                    Q[j] += λ 
                end
            end
            Q[i] = μ 
        else # Recovery event
            for j in neighbors(G, i)
                if x[j] == 0
                    Q[j] -= λ 
                end
            end
        end
        stats(t, x, i) 
    end
    return stats.vI, stats.vR 
end