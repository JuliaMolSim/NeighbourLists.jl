

module NeighbourLists

using StaticArrays

const SVec{T} = SVector{3, T}
const SMat{T} = SMatrix{3, 3, T, 9}

export CellList, npairs, nsites


# this contains the cell-list data structures and assembly
include("cell_list.jl")

# this contains the different iterators over sites, bonds, etc
include("iterators.jl")

# alternative assembly protocol more akin to mapreduce
include("mapreduce.jl")


end # module
