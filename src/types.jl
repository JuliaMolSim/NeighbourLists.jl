
using StaticArrays

const SVec{T} = SVector{3, T}
const SMat{T} = SMatrix{3, 3, T, 9}

export PairList, SortedCellList


"""
`PairList` stores a neighbourlist as a list of pairs.

### Constructors
```julia
PairList(X, cutoff, cell, pbc; backend=CPU())
PairList(sys::AbstractSystem, cutoff; backend=CPU())
```

### Fields
- `X` : positions as `Vector{SVec{T}}` (or GPU array)
- `C` : 3×3 cell matrix
- `cutoff` : cutoff distance
- `i, j` : pair indices, `(i[n], j[n])` is a neighbour pair
- `S` : cell shift vectors for periodic images
- `first` : CSR-like offsets, `first[m]:first[m+1]-1` are pairs for atom m
"""
struct PairList{T <: Real, TI <: Integer, AT <: AbstractVector}
   X::AT                        # Vector{SVec{T}} or CuVector{SVec{T}}
   C::SMat{T}
   cutoff::T
   i::AbstractVector{TI}
   j::AbstractVector{TI}
   S::AbstractVector{SVec{TI}}
   first::AbstractVector{TI}
end


"""
`SortedCellList` : GPU-friendly cell list using sort-based construction.

Atoms are sorted by their cell ID, enabling coalesced memory access.
`cell_offsets[c]:cell_offsets[c+1]-1` gives indices of atoms in cell `c`.

### Fields
- `X` : positions sorted by cell ID
- `perm` : permutation mapping sorted index → original index
- `cell_id` : cell ID for each atom (sorted)
- `cell_offsets` : CSR-like offsets for atoms in each cell
- `cell` : 3×3 cell matrix
- `inv_cell` : inverse of cell matrix
- `pbc` : periodic boundary conditions
- `cutoff` : cutoff distance
- `ncells` : number of cells in each dimension
- `ncells_total` : total number of cells
"""
struct SortedCellList{T <: Real, TI <: Integer, AT <: AbstractVector}
    X::AT                           # Positions sorted by cell (Vector{SVec{T}} or GPU)
    X_orig::AT                      # Original unsorted positions
    perm::AbstractVector{TI}        # sorted index → original index
    cell_id::AbstractVector{TI}     # Cell ID for each atom (sorted order)
    cell_offsets::AbstractVector{TI} # cell_offsets[c]:cell_offsets[c+1]-1 = atoms in cell c
    cell::SMat{T}
    inv_cell::SMat{T}
    pbc::SVec{Bool}
    cutoff::T
    ncells::SVec{TI}                # Number of cells per dimension
    ncells_total::TI                # Total number of cells
end


"""
`CellList` : (Legacy) store atoms in cells / bins using linked lists.
Mostly used internally to construct PairLists. Not GPU-compatible.
"""
struct CellList{T <: Real, TI <: Integer}
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
