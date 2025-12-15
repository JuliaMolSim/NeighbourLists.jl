# Tests for AtomsBase integration via the unified API
# neighbour_list(), neighbours(), build_cell_list()

using AtomsBase
using AtomsBuilder
using NeighbourLists
using NeighbourLists: neighbour_list, neighbours, num_neighbours, build_cell_list,
                      for_each_neighbour, SortedCellList
using Test
using Unitful

@testset "AtomsBase Integration" begin

    @testset "neighbour_list() - Periodic Bulk" begin
        sys = bulk(:Cu, cubic=true) * (4, 2, 3)
        cutoff = 3.5u"Å"

        # Unified API entry point
        nlist = neighbour_list(sys, cutoff)
        @test nlist isa PairList
        @test npairs(nlist) > 0
        @test nsites(nlist) == length(sys)

        # neighbours() function
        j, R = neighbours(nlist, 1)
        @test length(j) == 12  # FCC Cu has 12 nearest neighbours

        # num_neighbours() function
        @test num_neighbours(nlist, 1) == 12

        # All atoms in bulk Cu should have same coordination
        for i in 1:nsites(nlist)
            @test num_neighbours(nlist, i) == 12
        end
    end

    @testset "neighbour_list() - Lazy Mode" begin
        sys = bulk(:Cu, cubic=true) * (3, 3, 3)
        cutoff = 3.5u"Å"

        # Lazy mode returns SortedCellList
        clist = neighbour_list(sys, cutoff; lazy=true)
        @test clist isa SortedCellList
        @test nsites(clist) == length(sys)

        # neighbours() works with SortedCellList
        j, R, S = neighbours(clist, 1)
        @test length(j) == 12

        # for_each_neighbour works
        @test count_neighbours(clist, 1) == 12

        # Lazy and materialized give same pair count
        nlist = neighbour_list(sys, cutoff)
        @test npairs(nlist) == count_lazy_pairs(clist)
    end

    @testset "build_cell_list() - Direct" begin
        sys = bulk(:Cu, cubic=true) * (2, 2, 2)
        cutoff = 3.5u"Å"

        # Direct build_cell_list call
        clist = build_cell_list(sys, cutoff)
        @test clist isa SortedCellList
        @test nsites(clist) == length(sys)

        # Can iterate and get neighbours
        j, R, S = neighbours(clist, 1)
        @test length(j) == 12
    end

    @testset "Isolated System" begin
        isys = isolated_system([
            :H => [0., 0., 0.]u"Å",
            :H => [0., 0., 2.]u"Å",
            :H => [0., 10., 0.]u"Å"
        ])

        nlist = neighbour_list(isys, 5.0u"Å")

        # neighbours() function
        j, R = neighbours(nlist, 1)
        @test length(j) == 1
        @test j[1] == 2
        @test all(R[1] .≈ [0., 0., 2.])

        j, R = neighbours(nlist, 2)
        @test length(j) == 1
        @test j[1] == 1
        @test all(R[1] .≈ [0., 0., -2.])

        # Atom 3 is far from others
        j, R = neighbours(nlist, 3)
        @test length(j) == 0

        # num_neighbours()
        @test num_neighbours(nlist, 1) == 1
        @test num_neighbours(nlist, 2) == 1
        @test num_neighbours(nlist, 3) == 0

        # Lazy mode with isolated system
        clist = neighbour_list(isys, 5.0u"Å"; lazy=true)
        j_lazy, _, _ = neighbours(clist, 1)
        @test length(j_lazy) == 1
    end

    @testset "Different Cutoffs" begin
        sys = bulk(:Cu, cubic=true) * (3, 3, 3)

        # First nearest neighbours only (~2.55 Å for Cu)
        nlist_1nn = neighbour_list(sys, 2.6u"Å")
        @test num_neighbours(nlist_1nn, 1) == 12

        # Include second nearest neighbours (~3.6 Å)
        nlist_2nn = neighbour_list(sys, 3.7u"Å")
        @test num_neighbours(nlist_2nn, 1) > 12

        # Very small cutoff - no neighbours
        nlist_small = neighbour_list(sys, 1.0u"Å")
        @test npairs(nlist_small) == 0
    end

    @testset "Unit Handling" begin
        sys = bulk(:Cu, cubic=true) * (2, 2, 2)

        # Different unit specifications
        nlist_A = neighbour_list(sys, 3.5u"Å")
        nlist_nm = neighbour_list(sys, 0.35u"nm")
        nlist_pm = neighbour_list(sys, 350.0u"pm")

        # All should give same result
        @test npairs(nlist_A) == npairs(nlist_nm)
        @test npairs(nlist_A) == npairs(nlist_pm)
    end

    @testset "Error Handling" begin
        # 2D systems should throw error
        isys2d = isolated_system([
            :H => [0., 0.]u"Å",
            :H => [0., 2.]u"Å",
            :H => [10., 0.]u"Å"
        ])
        @test_throws ErrorException neighbour_list(isys2d, 5.0u"Å")
        @test_throws ErrorException build_cell_list(isys2d, 5.0u"Å")
    end

    @testset "Legacy PairList Constructor" begin
        # Legacy API should still work
        sys = bulk(:Cu, cubic=true) * (4, 2, 3)
        nlist = PairList(sys, 3.5u"Å")

        j, R = neighbours(nlist, 1)
        @test length(j) == 12
    end
end
