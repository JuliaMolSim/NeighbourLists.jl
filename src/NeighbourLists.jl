
__precompile__()
module NeighbourLists

using KernelAbstractions
using KernelAbstractions: @kernel, @index, @Const, synchronize, get_backend
using AcceleratedKernels

# Re-export backend types for user convenience
# I don't think this should be exported, it is too generic a name. 
# export CPU

# Default backend
default_backend() = CPU()

# Get backend from array (auto-detect GPU arrays)
get_array_backend(X::AbstractArray) = get_backend(X)
get_array_backend(::Any) = CPU()

# Legacy threading control (for old cell_list code)
const MAX_THREADS = [1]
set_maxthreads!(n) = (MAX_THREADS[1] = n)

include("types.jl")

# Adapt.adapt_structure rules for SortedCellList / PairList so external
# packages can pass them as KA kernel arguments (fixes #37).
include("adapt.jl")

# this contains the cell-list data structures and assembly
include("cell_list.jl")

# GPU kernels (portable via KernelAbstractions)
include("gpu_kernels.jl")

# this contains the different iterators over sites, bonds, etc
include("iterators.jl")

# AtomsBase interface is provided via extension: ext/NeighbourListsAtomsBaseExt.jl
# CUDA support is provided via extension: ext/NeighbourListsCUDAExt.jl

end # module
