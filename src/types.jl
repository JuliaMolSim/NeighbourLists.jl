
using StaticArrays

const SVec{T} = SVector{3, T}
const SMat{T} = SMatrix{3, 3, T, 9}

export PairList


"""
`PairList` stores a neighbourlist as a list of pairs
```
PairList(X, cutoff, cell, pbc)
PairList(nlist::CellList)
```
where

* `X` : positions, either as 3 x N matrix or as `Vector{SVec}`
* `cutoff` : positive real value
* `cell` : 3 x 3 matrix, with rows denoting the cell vectors
* `pbc` : 3-tuple or array, storing periodicity information

### Kw-args:

* `int_type` : default is `Int`
* `store_first` : whether to store the array of first indices, default `true`
* `sorted` : whether to sort the `j` vector, default `false`

### PairList fields

`X, cutoff, i, j, r, R, first`, where

`(i[n], j[n])` denotes the indices of a neighbour pair, `r[n]` the distance
between those atoms, `R[n]` the vectorial distance, note this is identical to
`X[i[n]]-X[j[n]]` without periodic b.c.s, but with periodic boundary conditions
it is different. `first[m]` contains the index to the first `(i[n], j[n])` for
which `i[n] == first[m]`, i.e., `(j, first)` essentially defines a compressed
column storage of the adjacancy matrix.
"""
struct PairList{T <: AbstractFloat, TI <: Integer}
   X::Vector{SVec{T}}
   C::SMat{T}
   cutoff::T
   i::Vector{TI}
   j::Vector{TI}
   S::Vector{SVec{TI}}
   first::Vector{TI}
end


"""
`CellList` : store atoms in cells / bins. Mostly used internally
to construct PairLists.
"""
struct CellList{T <: AbstractFloat, TI <: Integer}
   X::Vector{SVec{T}}
   cell::SMat{T}
   inv_cell::SMat{T}
   pbc::SVec{Bool}
   cutoff::T
   seed::Vector{TI}
   last::Vector{TI}
   next::Vector{TI}
   nats::Vector{TI}
end
