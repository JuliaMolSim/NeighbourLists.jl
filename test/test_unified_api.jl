# Tests for the unified high-level API (neighbour_list, neighbours, num_neighbours)
# These functions are the recommended entry points for users

using NeighbourLists: neighbour_list, neighbours, num_neighbours, 
                      max_neighbours, maxneigs, CPU, neigss 

@testset "Unified API" begin

    @testset "neighbour_list() - Materialized" begin
        X, C, L = rand_config(100)
        cutoff = L / 3

        # Basic usage - should return PairList
        nlist = neighbour_list(X, cutoff, C, FULL_PBC)
        @test nlist isa PairList
        @test npairs(nlist) > 0
        @test nsites(nlist) == 100

        # Verify against direct build_cell_list + materialize
        clist = build_cell_list(X, cutoff, C, FULL_PBC; backend=CPU())
        nlist_direct = materialize_pairlist(clist)
        @test npairs(nlist) == npairs(nlist_direct)
        @test compare_pairlists(nlist, nlist_direct)
    end

    @testset "neighbour_list() - Lazy Mode" begin
        X, C, L = rand_config(100)
        cutoff = L / 3

        # Lazy mode should return SortedCellList
        clist = neighbour_list(X, cutoff, C, FULL_PBC; lazy=true)
        @test clist isa SortedCellList
        @test nsites(clist) == 100

        # Lazy and materialized should give same results
        nlist = neighbour_list(X, cutoff, C, FULL_PBC; lazy=false)

        # Count pairs via lazy iteration - should match materialized count
        @test count_lazy_pairs(clist) == npairs(nlist)
    end

    @testset "neighbour_list() - Integer Types" begin
        X, C, L = rand_config(50)
        cutoff = L / 3

        for TI in [Int32, Int64]
            nlist = neighbour_list(X, cutoff, C, FULL_PBC; int_type=TI)
            @test eltype(nlist.i) == TI
            @test eltype(nlist.j) == TI
        end
    end

    @testset "neighbour_list() - All PBC Combinations" begin
        X, C, L = rand_config(50)
        cutoff = L / 4

        for pbc in all_pbc_cases()
            nlist = neighbour_list(X, cutoff, C, pbc)
            @test nlist isa PairList

            # Verify against legacy implementation
            nlist_legacy = PairList(X, cutoff, C, pbc; int_type=Int32)
            @test npairs(nlist) == npairs(nlist_legacy)
        end
    end

    @testset "neighbour_list() - Edge Cases" begin
        C = cubic_cell(10.0)

        # Empty system - single atom
        nlist = neighbour_list([SVec(5.0, 5.0, 5.0)], 3.0, C, FULL_PBC)
        @test npairs(nlist) == 0

        # Two atoms within cutoff
        X = [SVec(5.0, 5.0, 5.0), SVec(5.0, 5.0, 6.0)]
        nlist = neighbour_list(X, 3.0, C, FULL_PBC)
        @test npairs(nlist) == 2  # Both directions

        # Two atoms outside cutoff (no PBC)
        X = [SVec(1.0, 1.0, 1.0), SVec(8.0, 8.0, 8.0)]
        nlist = neighbour_list(X, 3.0, C, NO_PBC)
        @test npairs(nlist) == 0
    end
end

@testset "neighbours() Function" begin

    @testset "neighbours() with PairList" begin
        X, C, L = rand_config(50)
        nlist = neighbour_list(X, L/3, C, FULL_PBC)

        for i in 1:5
            # neighbours() should work as alias for neigs()
            j_neigh, R_neigh, S_neigh = neighbours(nlist, i)
            j_neigs, R_neigs = neigs(nlist, i)
            j1_neigs, R1_neigs, S1_neigs = neigss(nlist, i)

            @test j_neigh == j_neigs == j1_neigs
            @test S_neigh == S1_neigs
            @test R_neigh == R_neigs == R1_neigs
        end
    end

    @testset "neighbours() Consistency Between PairList and SortedCellList" begin
        X, C, L = rand_config(50)
        nlist = neighbour_list(X, L/3, C, FULL_PBC)
        clist = neighbour_list(X, L/3, C, FULL_PBC; lazy=true)

        # Should have same neighbours for all atoms (order may differ)
        @test verify_neighbours_consistency(nlist, clist)
    end
end

@testset "num_neighbours() Function" begin

    @testset "num_neighbours() with PairList" begin
        X, C, L = rand_config(50)
        nlist = neighbour_list(X, L/3, C, FULL_PBC)

        # Count manually and verify
        expected_counts = count_manual_neighbours(nlist)
        for i in 1:50
            @test num_neighbours(nlist, i) == expected_counts[i]
        end
    end

    @testset "num_neighbours() with SortedCellList" begin
        X, C, L = rand_config(50)
        clist = neighbour_list(X, L/3, C, FULL_PBC; lazy=true)

        for i in 1:50
            # num_neighbours should match count_neighbours
            @test num_neighbours(clist, i) == count_neighbours(clist, i)
        end
    end

    @testset "num_neighbours() Consistency" begin
        X, C, L = rand_config(50)
        nlist = neighbour_list(X, L/3, C, FULL_PBC)
        clist = neighbour_list(X, L/3, C, FULL_PBC; lazy=true)

        for i in 1:50
            @test num_neighbours(nlist, i) == num_neighbours(clist, i)
        end
    end
end

@testset "max_neighbours() Function" begin

    @testset "max_neighbours() Basic" begin
        X, C, L = rand_config(100)
        nlist = neighbour_list(X, L/3, C, FULL_PBC)

        max_n = max_neighbours(nlist)
        @test max_n >= 0
        @test max_n == maxneigs(nlist)  # Should be same as maxneigs alias

        # Verify it's actually the maximum
        for i in 1:nsites(nlist)
            @test num_neighbours(nlist, i) <= max_n
        end
    end

    @testset "max_neighbours() Edge Cases" begin
        C = cubic_cell(10.0)

        # Single atom - max neighbours is 0
        nlist = neighbour_list([SVec(5.0, 5.0, 5.0)], 3.0, C, FULL_PBC)
        @test max_neighbours(nlist) == 0

        # Two atoms - each has 1 neighbour
        X = [SVec(5.0, 5.0, 5.0), SVec(5.0, 5.0, 6.0)]
        nlist = neighbour_list(X, 3.0, C, FULL_PBC)
        @test max_neighbours(nlist) == 1
    end
end

@testset "GPU - Unified API" begin
    if gpu_available()
        backend_name = gpu_backend()
        # use Float32 on Metal (no Float64 support), Float64 elsewhere
        Tgpu = first(gpu_supported_eltypes())
        @info "Testing unified API with $backend_name ($Tgpu)"

        @testset "neighbour_list() with GPU ($backend_name)" begin
            X, C, L = rand_config(100; T=Tgpu)
            X_gpu = to_gpu_array(X)

            # Auto-detect GPU backend from array type
            nlist_gpu = neighbour_list(X_gpu, L/3, C, FULL_PBC)
            @test nlist_gpu isa PairList
            @test nlist_gpu.i isa gpu_array_type()
            @test npairs(nlist_gpu) > 0

            # Compare with CPU
            nlist_cpu = neighbour_list(X, L/3, C, FULL_PBC)
            @test npairs(nlist_cpu) == npairs(nlist_gpu)
            @test compare_cpu_gpu_full(nlist_cpu, nlist_gpu)
        end

        @testset "neighbour_list() Lazy Mode GPU ($backend_name)" begin
            X, C, L = rand_config(100; T=Tgpu)
            X_gpu = to_gpu_array(X)

            clist_gpu = neighbour_list(X_gpu, L/3, C, FULL_PBC; lazy=true)
            @test clist_gpu isa SortedCellList
            @test clist_gpu.X_orig isa gpu_array_type()
        end
    else
        @testset "GPU Unified API (Skipped)" begin
            @test_skip "No GPU backend available"
        end
    end
end
