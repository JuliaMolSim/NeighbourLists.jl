using AtomsBase
using StaticArrays: SVector
using Unitful


function PairList(ab::AtomsBase.AbstractSystem, cutoff::Unitful.Length; length_unit=unit(cutoff))
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
    
    nlist = PairList(r, ustrip(length_unit, cutoff), cell_matrix, pbc; int_type=Int)
    return nlist
end