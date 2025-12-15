# Tests for sort-based cell list implementation
# Uses shared utilities from test_utils.jl (included by runtests.jl)

using NeighbourLists: SortedCellList, map_sites!, map_pairs!, map_pairs_d!,
                      count_neighbours, get_neighbours

@testset "Sort-based Cell List" begin

    @testset "Basic Construction" begin
        X, C, L = rand_config(100)
        clist = build_cell_list(X, L/3, C, FULL_PBC; backend=CPU())

        @test clist isa SortedCellList
        @test nsites(clist) == 100
        @test npairs(materialize_pairlist(clist)) > 0
    end

    @testset "Legacy vs Sort-based" begin
        for _ in 1:5
            test_legacy_vs_sortbased(N=rand(50:200))
        end
    end

    @testset "All PBC Combinations" begin
        test_all_pbc_legacy_vs_sortbased()
    end

    @testset "Non-cubic Cells" begin
        test_triclinic_cell() do X, C, cutoff, pbc
            nlist_legacy = PairList(X, cutoff, C, pbc; int_type=Int32)
            nlist_sort = materialize_pairlist(build_cell_list(X, cutoff, C, pbc; backend=CPU()))
            @test npairs(nlist_legacy) == npairs(nlist_sort)
            @test compare_pairlists(nlist_legacy, nlist_sort)
        end

        # Elongated cell
        C2 = SMat(diagm([5.0, 5.0, 20.0]))
        X2 = [C2' * SVec(rand(), rand(), rand()) for _ in 1:80]
        nlist1 = PairList(X2, 3.0, C2, FULL_PBC; int_type=Int32)
        nlist2 = materialize_pairlist(build_cell_list(X2, 3.0, C2, FULL_PBC; backend=CPU()))
        @test compare_pairlists(nlist1, nlist2)
    end

    @testset "Edge Cases" begin
        test_edge_cases() do X, cutoff, C, pbc
            materialize_pairlist(build_cell_list(X, cutoff, C, pbc; backend=CPU()))
        end
    end

    @testset "Large Cutoff" begin
        X, C, L = rand_config(30)
        cutoff = L * 0.6  # > L/2
        nlist1 = PairList(X, cutoff, C, FULL_PBC; int_type=Int32)
        nlist2 = materialize_pairlist(build_cell_list(X, cutoff, C, FULL_PBC; backend=CPU()))
        @test compare_pairlists(nlist1, nlist2)
    end

    @testset "Integer Types" begin
        X, C, L = rand_config(100)
        for TI in [Int32, Int64]
            nlist = materialize_pairlist(build_cell_list(X, L/3, C, FULL_PBC; backend=CPU(), int_type=TI))
            @test eltype(nlist.i) == TI
        end
    end

    @testset "Large Systems" begin
        test_sizes(sizes=[500, 1000, 2000]) do X, C, L, cutoff
            nlist1 = PairList(X, cutoff, C, FULL_PBC; int_type=Int32)
            nlist2 = materialize_pairlist(build_cell_list(X, cutoff, C, FULL_PBC; backend=CPU()))
            @test npairs(nlist1) == npairs(nlist2)
        end
    end
end

@testset "Lazy Iteration" begin
    X, C, L = rand_config(100)
    clist = build_cell_list(X, L/3, C, FULL_PBC; backend=CPU())
    nlist = materialize_pairlist(clist)

    # Count via lazy iteration
    lazy_count = sum(1:nsites(clist)) do i
        count = 0
        for_each_neighbour(clist, i) do j, R, S
            count += 1
        end
        count
    end
    @test lazy_count == npairs(nlist)

    # Collect pairs via lazy iteration
    lazy_pairs = Tuple{Int, Int, NTuple{3,Int}}[]
    for i in 1:nsites(clist)
        for_each_neighbour(clist, i) do j, R, S
            push!(lazy_pairs, (i, j, Tuple(S)))
        end
    end
    @test sort(lazy_pairs) == get_sorted_pairs(nlist)
end

@testset "get_neighbours / count_neighbours" begin
    X, C, L = rand_config(50)
    clist = build_cell_list(X, L/3, C, FULL_PBC; backend=CPU())
    nlist = materialize_pairlist(clist)

    # Build expected counts from PairList
    expected = zeros(Int, 50)
    for idx in 1:npairs(nlist)
        expected[nlist.i[idx]] += 1
    end

    for i in 1:50
        @test count_neighbours(clist, i) == expected[i]
        js, Rs, Ss = get_neighbours(clist, i)
        @test length(js) == expected[i]
    end
end

@testset "MapReduce Functions" begin
    X, C, L = rand_config(100)
    clist = build_cell_list(X, L/3, C, FULL_PBC; backend=CPU())
    nlist = materialize_pairlist(clist)

    @testset "map_sites!" begin
        out = zeros(Float64, 100)
        map_sites!(Rs -> sum(norm.(Rs)), out, clist)
        expected = zeros(Float64, 100)
        for i in 1:100
            for_each_neighbour(clist, i) do j, R, S
                expected[i] += norm(R)
            end
        end
        @test out ≈ expected atol=1e-10
    end

    @testset "map_pairs! symmetric" begin
        out = zeros(Float64, 100)
        map_pairs!((i, j, R) -> 1.0 / norm(R), out, nlist; symmetric=true)
        # Verify it ran without error and produced output
        @test sum(out) > 0
    end

    @testset "map_pairs_d! anti-symmetric" begin
        out = zeros(SVec{Float64}, 100)
        map_pairs_d!((i, j, R) -> R / norm(R), out, nlist)
        # Anti-symmetric: total should be ~zero
        @test norm(sum(out)) < 1e-10
    end
end
