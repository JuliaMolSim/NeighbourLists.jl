

module NeighbourList

using StaticArrays

const SVec{T} = SVector{3, T}
const SMat{T} = SMatrix{3, 3, T, 9}

export CellList, npairs, nsites


# this contains the neighourlist assemby routine
include("core.jl")

# ==================== GateWay Routines  ================

"""
Basic `NeighbourList` type. Typically, this is constructed using
"""
struct CellList{T <: AbstractFloat, TI <: Integer}
   X::Vector{SVec{T}}
   cutoff::T
   i::Vector{TI}
   j::Vector{TI}
   r::Vector{T}
   R::Vector{SVec{T}}
   S::Vector{SVec{TI}}
   first::Vector{TI}
end

# TODO: consider redefining the cell in terms of the
#       extent of X

function CellList{T}(X::Vector{SVec{T}}, cutoff::AbstractFloat,
                     cell::AbstractMatrix, pbc;
                     int_type::Type = Int)
   i, j, r, R, S = _neighbour_list_(SMat{T}(cell...), SVec{Bool}(pbc...), X,
                                    cutoff, zero(int_type))
   return CellList(X, cutoff, i, j, r, R, S, get_first(i, length(X)))
end

CellList{T}(X::Matrix{T}, args...; kwargs...) =
      CellList(reinterpret(SVec{T}, X, (size(X,2),)), args...; varargs...)

npairs(nlist::CellList) = length(nlist.i)
nsites(nlist::CellList) = length(nlist.first) - 1

# this contains the different iterators over sites, bonds, etc
include("iterators.jl")
include("reduce.jl")

end # module
