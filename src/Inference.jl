dist(A, B) = norm(A .- B)

function likelihood(G, λ, μ, T, M, β, Ostar)
    L = 0.0
    vT = 1:(T/20):T 
    x = zeros(Int, size(G, 1)) # Pre-allocated to avoid global variable type instability
    
    for m in 1:M
        stats = DeltaRecorder(vT) # Reset stats cleanly for each sample
        temp_I, temp_R = sample!(x, G, λ, μ, T; stats=stats)
        
        O_I = cumsum(temp_I)
        O_R = cumsum(temp_R)
        
        dist_I = dist(O_I, Ostar[1])
        # dist_R = dist(O_R, Ostar[2]) 
        total_dist = dist_I # + dist_R
        
        L += exp(-β * total_dist)
    end
    return L / M
end

function maxlikelihood(λ_range, μ_range, G, T, M, Ostar, β)
    Lmax = zeros(length(λ_range), length(μ_range))
    for i in 1:length(λ_range) 
        for r in 1:length(μ_range)
            Lmax[i, r] = likelihood(G, λ_range[i], μ_range[r], T, M, β, Ostar)
        end
    end
    return Lmax
end

function log_likelihood(G, λ, μ, T, M, β, Ostar)
    log_weights = fill(0.0, M)
    vT = 1:(T/20):T
    x = zeros(Int, size(G, 1))
    
    for m in 1:M
        stats = DeltaRecorder(vT) 
        temp_I, temp_R = sample!(x, G, λ, μ, T; stats=stats)
        
        O_I = cumsum(temp_I)
        O_R = cumsum(temp_R)
        
        dist_I = dist(O_I, Ostar[1])
        # dist_R = dist(O_R, Ostar[2])
        total_dist = dist_I # + dist_R
        
        log_weights[m] = -β * total_dist
    end

    max_val = maximum(log_weights)
    L = -log(M) + max_val + log(sum(exp.(log_weights .- max_val)))
    
    return L
end

function maxloglikelihood(λ_range, μ_range, G, T, M, Ostar, β)
    Lmax = zeros(length(λ_range), length(μ_range))
    for i in 1:length(λ_range)
        for r in 1:length(μ_range)
            Lmax[i, r] = log_likelihood(G, λ_range[i], μ_range[r], T, M, β, Ostar)
        end
    end
    return Lmax
end