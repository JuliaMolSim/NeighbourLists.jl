
__precompile__()
module NeighbourLists

const MAX_THREADS = [1]   # temporarily hard-coded single thread!
set_maxthreads!(n) = (NeighbourLists.MAX_THREADS[1] = n)


include("types.jl")

# this contains the cell-list data structures and assembly
include("cell_list.jl")

# this contains the different iterators over sites, bonds, etc
include("iterators.jl")

# alternative assembly protocol more akin to mapreduce
include("mapreduce.jl")


end # module
