using AtomsBase
using StaticArrays: SVector
using Unitful


function PairList(ab::AtomsBase.AbstractSystem, cutoff::Unitful.Length; length_unit=unit(cutoff))
    cell = ustrip.(length_unit, hcat( bounding_box(ab)... )' )
    pbc = periodicity(ab)
    r = map( 1:length(ab) ) do i
        # Need to have SVector here for PairList to work
        # if position does not give SVector
        SVector( ustrip.(length_unit, position(ab,i))...)
    end
    nlist = PairList(r, ustrip(length_unit, cutoff), cell, pbc; int_type=Int)
    return nlist
end