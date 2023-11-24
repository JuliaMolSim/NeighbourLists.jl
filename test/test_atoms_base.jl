using ASEconvert
using NeighbourLists
using Test
using Unitful

@testset "AtomsBase PairList" begin
    cu = ase.build.bulk("Cu") * pytuple((4, 2, 3))
    sys = pyconvert(AbstractSystem, cu)

    nlist = PairList(sys, 3.5u"Å")

    id, r = neigs(nlist,1)

    @test length(id) == 12
end




