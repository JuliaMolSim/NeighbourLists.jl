module NeighbourListsCUDAExt

# This extension exists to ensure CUDA.jl is properly loaded when available.
# The actual GPU kernels are portable via KernelAbstractions.jl + AcceleratedKernels.jl
# and work across all GPU backends (CUDA, ROCm, Metal, oneAPI).
#
# Having this as a weak dependency means users only need to `using CUDA` to get
# GPU acceleration - no additional NeighbourLists-specific setup required.

using NeighbourLists
using CUDA

end # module
