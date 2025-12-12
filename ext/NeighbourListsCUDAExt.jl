module NeighbourListsCUDAExt

using NeighbourLists
using CUDA

# CUDA extension for NeighbourLists.jl
# The main GPU code is already portable via KernelAbstractions.jl
# This extension ensures CUDA is properly loaded when available

# Re-export CuArray for convenience
# Users can do: using NeighbourLists, CUDA; X_gpu = CuArray(X)

end # module
