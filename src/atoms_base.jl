using AtomsBase
using StaticArrays: SVector
using Unitful


"""
    PairList(system::AbstractSystem, cutoff::Unitful.Length; kwargs...)

Construct a PairList from an AtomsBase system.

# Keyword arguments
- `length_unit`: Unit for positions (default: unit of cutoff)
- `backend`: KernelAbstractions backend (default: CPU())
- `int_type`: Integer type for indices (default: Int)
"""
function PairList(ab::AtomsBase.AbstractSystem, cutoff::Unitful.Length;
                  length_unit=unit(cutoff),
                  backend=default_backend(),
                  int_type::Type=Int)
    r = map( 1:length(ab) ) do i
        # Need to have SVector here for PairList to work
        # if position does not give SVector
        SVector( ustrip.(length_unit, position(ab,i))...)
    end
    tmp_cell = cell(ab)
    if tmp_cell isa AtomsBase.IsolatedCell{3, <:Any}
        # this only works in 3D
        max_x = maximum(q->q[1], r)
        max_y = maximum(q->q[2], r)
        max_z = maximum(q->q[3], r)
        min_x = minimum(q->q[1], r)
        min_y = minimum(q->q[2], r)
        min_z = minimum(q->q[3], r)
        cell_matrix = SMatrix{3,3, typeof(max_x), 9}(
            # add +1 for cell size
            max_x - min_x + 1, 0., 0.,
            0., max_y - min_y + 1, 0.,
            0., 0., max_z - min_z + 1,
        )
    elseif tmp_cell isa AtomsBase.IsolatedCell
        D = n_dimensions(ab)
        throw( error("NeighbourLists.jl does not support $D-dimensional isolated AtomsBase systems yet.") )
    else
        cell_matrix = ustrip.(length_unit, hcat( cell_vectors(ab)... )' )
    end
    pbc = periodicity(ab)

    # Use the new sort-based cell list construction
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
function build_cell_list(ab::AtomsBase.AbstractSystem, cutoff::Unitful.Length;
                         length_unit=unit(cutoff),
                         backend=default_backend(),
                         int_type::Type=Int32)
    r = map( 1:length(ab) ) do i
        SVector( ustrip.(length_unit, position(ab,i))...)
    end
    tmp_cell = cell(ab)
    if tmp_cell isa AtomsBase.IsolatedCell{3, <:Any}
        max_x = maximum(q->q[1], r)
        max_y = maximum(q->q[2], r)
        max_z = maximum(q->q[3], r)
        min_x = minimum(q->q[1], r)
        min_y = minimum(q->q[2], r)
        min_z = minimum(q->q[3], r)
        cell_matrix = SMatrix{3,3, typeof(max_x), 9}(
            max_x - min_x + 1, 0., 0.,
            0., max_y - min_y + 1, 0.,
            0., 0., max_z - min_z + 1,
        )
    elseif tmp_cell isa AtomsBase.IsolatedCell
        D = n_dimensions(ab)
        throw( error("NeighbourLists.jl does not support $D-dimensional isolated AtomsBase systems yet.") )
    else
        cell_matrix = ustrip.(length_unit, hcat( cell_vectors(ab)... )' )
    end
    pbc = periodicity(ab)

    return build_cell_list(r, ustrip(length_unit, cutoff), cell_matrix, pbc;
                           int_type=int_type, backend=backend)
end