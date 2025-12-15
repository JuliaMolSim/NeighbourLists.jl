module NeighbourListsAtomsBaseExt

using NeighbourLists
using NeighbourLists: SMat, SVec, default_backend, build_cell_list, materialize_pairlist,
                      neighbour_list, neighbours, bool_to_val
using AtomsBase
using AtomsBase: AbstractSystem, IsolatedCell, cell, cell_vectors, position, periodicity, n_dimensions
using StaticArrays: SVector, SMatrix
using Unitful

"""
    _get_cell_matrix(ab, positions, length_unit)

Extract cell matrix from AtomsBase system, handling IsolatedCell specially.
For isolated systems, constructs a bounding box from positions.
"""
function _get_cell_matrix(ab::AbstractSystem, positions, length_unit)
    tmp_cell = cell(ab)
    if tmp_cell isa IsolatedCell{3, <:Any}
        # For isolated 3D systems, construct bounding box
        max_x = maximum(q->q[1], positions)
        max_y = maximum(q->q[2], positions)
        max_z = maximum(q->q[3], positions)
        min_x = minimum(q->q[1], positions)
        min_y = minimum(q->q[2], positions)
        min_z = minimum(q->q[3], positions)
        return SMatrix{3,3, typeof(max_x), 9}(
            max_x - min_x + 1, 0., 0.,
            0., max_y - min_y + 1, 0.,
            0., 0., max_z - min_z + 1,
        )
    elseif tmp_cell isa IsolatedCell
        D = n_dimensions(ab)
        throw(error("NeighbourLists.jl does not support $D-dimensional isolated AtomsBase systems yet."))
    else
        return ustrip.(length_unit, hcat(cell_vectors(ab)...)')
    end
end

"""
    PairList(system::AbstractSystem, cutoff::Unitful.Length; kwargs...)

Construct a PairList from an AtomsBase system.

# Keyword arguments
- `length_unit`: Unit for positions (default: unit of cutoff)
- `backend`: KernelAbstractions backend (default: CPU())
- `int_type`: Integer type for indices (default: Int32)
"""
function NeighbourLists.PairList(ab::AbstractSystem, cutoff::Unitful.Length;
                                  length_unit=unit(cutoff),
                                  backend=default_backend(),
                                  int_type::Type=Int32)
    r = map(1:length(ab)) do i
        SVector(ustrip.(length_unit, position(ab, i))...)
    end
    cell_matrix = _get_cell_matrix(ab, r, length_unit)
    pbc = periodicity(ab)

    clist = build_cell_list(r, ustrip(length_unit, cutoff), cell_matrix, pbc;
                            int_type=int_type, backend=backend)
    return materialize_pairlist(clist; backend=backend)
end

"""
    build_cell_list(system::AbstractSystem, cutoff::Unitful.Length; kwargs...)

Build a SortedCellList from an AtomsBase system (without materializing pairs).

# Keyword arguments
- `length_unit`: Unit for positions (default: unit of cutoff)
- `backend`: KernelAbstractions backend (default: CPU())
- `int_type`: Integer type for indices (default: Int32)
"""
function NeighbourLists.build_cell_list(ab::AbstractSystem, cutoff::Unitful.Length;
                                         length_unit=unit(cutoff),
                                         backend=default_backend(),
                                         int_type::Type=Int32)
    r = map(1:length(ab)) do i
        SVector(ustrip.(length_unit, position(ab, i))...)
    end
    cell_matrix = _get_cell_matrix(ab, r, length_unit)
    pbc = periodicity(ab)

    return build_cell_list(r, ustrip(length_unit, cutoff), cell_matrix, pbc;
                           int_type=int_type, backend=backend)
end

"""
    neighbour_list(system::AbstractSystem, cutoff::Unitful.Length; kwargs...)

Unified entry point for building neighbour lists from AtomsBase systems.

Returns a `PairList` by default. If `lazy=true`, returns a `SortedCellList`
that can be iterated without materializing all pairs (memory efficient).

# Keyword arguments
- `length_unit`: Unit for positions (default: unit of cutoff)
- `backend`: KernelAbstractions backend (default: CPU())
- `lazy`: If true, return SortedCellList instead of materializing pairs
- `int_type`: Integer type for indices (default: Int32)

# Examples
```julia
using AtomsBuilder, Unitful

sys = bulk(:Cu, cubic=true) * (3,3,3)

# Materialized pair list (default)
nlist = neighbour_list(sys, 5.0u"Å")
for i in 1:length(sys)
    j, R, S = neighbours(nlist, i)
    # process neighbours...
end

# Lazy cell list (memory efficient)
clist = neighbour_list(sys, 5.0u"Å"; lazy=true)
```
"""
function NeighbourLists.neighbour_list(ab::AbstractSystem, cutoff::Unitful.Length;
                                        length_unit=unit(cutoff),
                                        backend=default_backend(),
                                        lazy::Union{Bool, Val}=false,
                                        int_type::Type=Int32)
    return _neighbour_list_ab(bool_to_val(lazy), ab, cutoff; length_unit=length_unit,
                              backend=backend, int_type=int_type)
end

function _neighbour_list_ab(lazy::Val{L}, ab::AbstractSystem, cutoff::Unitful.Length;
                            length_unit=unit(cutoff),
                            backend=default_backend(),
                            int_type::Type=Int32) where {L}
    clist = build_cell_list(ab, cutoff; length_unit=length_unit,
                            backend=backend, int_type=int_type)
    if L
        return clist
    else
        return materialize_pairlist(clist; backend=backend)
    end
end

end # module
