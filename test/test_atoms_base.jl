#using ASEconvert
using AtomsBase
using AtomsBuilder
using NeighbourLists
using Test
using Unitful

@testset "AtomsBase PairList" begin
    sys = bulk(:Cu, cubic=true) * (4,2,3)

    nlist = PairList(sys, 3.5u"Å")

    id, r = neigs(nlist,1)

    @test length(id) == 12

    @testset "AtomsBase Isolated System" begin
        isys = isolated_system([
            :H => [0., 0., 0.]u"Å",
            :H => [0., 0., 2.]u"Å",
            :H => [0., 10., 0.]u"Å"
        ])
        nlist = PairList(isys, 5.0u"Å")
        id, r = neigs(nlist,1)
        @test length(id) == 1
        @test id[1] == 2
        @test all( r[1] .≈ [0., 0., 2.] )
        id, r = neigs(nlist,2)
        @test length(id) == 1
        @test id[1] == 1
        @test all( r[1] .≈ [0., 0., -2.] )
        id, r = neigs(nlist,3)
        @test length(id) == 0

        isys2d = isolated_system([
            :H => [ 0., 0.]u"Å",
            :H => [ 0., 2.]u"Å",
            :H => [ 10., 0.]u"Å"
        ])
        @test_throws ErrorException PairList(isys2d, 5.0u"Å")
    end
end




