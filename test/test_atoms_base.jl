#using ASEconvert
using AtomsBuilder
using NeighbourLists
using Test
using Unitful

@testset "AtomsBase PairList" begin
    sys = bulk(:Cu, cubic=true) * (4,2,3)

    nlist = PairList(sys, 3.5u"Å")

    id, r = neigs(nlist,1)

    @test length(id) == 12
end




