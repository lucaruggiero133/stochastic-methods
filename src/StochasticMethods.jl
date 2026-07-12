module StochasticMethods

using SparseArrays
using LinearAlgebra
using ExponentialQueues

# Export the functions we want available to the user
export generate_sparse_graph, sample!, DeltaRecorder
export dist, likelihood, maxlikelihood, log_likelihood, maxloglikelihood

# Include the sub-modules
include("Epidemic.jl")
include("Inference.jl")

end