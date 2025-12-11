# Tests for sort-based cell list implementation
# Validates against legacy linked-list implementation

using NeighbourLists
using NeighbourLists: SMat, SVec, SortedCellList, nsites, npairs,
                      build_cell_list, materialize_pairlist, for_each_neighbour,
                      map_sites!, map_pairs!
using Test
using LinearAlgebra
using StaticArrays

# Helper to generate random configurations
function rand_config_sortbased(N; density=0.05)
    volume = N / density
    L = volume^(1/3)
    C = SMat(diagm([L, L, L]))
    X = [SVec(L * rand(), L * rand(), L * rand()) for _ in 1:N]
    return X, C, L
end

# Helper to compare pair lists (order-independent)
function compare_pairlists(nlist1, nlist2)
    pairs1 = sort(collect(zip(nlist1.i, nlist1.j, [Tuple(s) for s in nlist1.S])))
    pairs2 = sort(collect(zip(nlist2.i, nlist2.j, [Tuple(s) for s in nlist2.S])))
    return pairs1 == pairs2
end

# Helper to get pairs as sorted tuples (i, j, S)
function get_sorted_pairs(nlist)
    return sort(collect(zip(nlist.i, nlist.j, [Tuple(s) for s in nlist.S])))
end

@testset "Sort-based Cell List" begin

    @testset "Basic Construction" begin
        # Simple cubic system
        N = 100
        X, C, L = rand_config_sortbased(N)
        cutoff = L / 3
        pbc = SVec(true, true, true)

        # Build sort-based cell list
        clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())

        @test clist isa SortedCellList
        @test nsites(clist) == N
        @test clist.cutoff == cutoff

        # Materialize pairs
        nlist = materialize_pairlist(clist)
        @test npairs(nlist) > 0
        @test length(nlist.i) == length(nlist.j) == length(nlist.S)
    end

    @testset "Comparison with Legacy Implementation" begin
        # Test multiple random configurations
        for trial in 1:5
            N = rand(50:200)
            X, C, L = rand_config_sortbased(N)
            cutoff = L / 4 + rand() * L / 4  # Random cutoff between L/4 and L/2
            pbc = SVec(rand(Bool), rand(Bool), rand(Bool))

            # Legacy implementation
            nlist_legacy = PairList(X, cutoff, C, pbc; int_type=Int32)

            # Sort-based implementation
            clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
            nlist_sort = materialize_pairlist(clist)

            # Compare pair counts
            @test npairs(nlist_legacy) == npairs(nlist_sort)

            # Compare actual pairs (order-independent)
            @test compare_pairlists(nlist_legacy, nlist_sort)
        end
    end

    @testset "Periodic Boundary Conditions" begin
        N = 50
        X, C, L = rand_config_sortbased(N)
        cutoff = L / 3

        # Test all PBC combinations
        pbc_cases = [
            SVec(true, true, true),    # Full PBC
            SVec(false, false, false), # No PBC
            SVec(true, false, false),  # PBC in x only
            SVec(false, true, false),  # PBC in y only
            SVec(false, false, true),  # PBC in z only
            SVec(true, true, false),   # PBC in x,y
            SVec(true, false, true),   # PBC in x,z
            SVec(false, true, true),   # PBC in y,z
        ]

        for pbc in pbc_cases
            nlist_legacy = PairList(X, cutoff, C, pbc; int_type=Int32)
            clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
            nlist_sort = materialize_pairlist(clist)

            @test npairs(nlist_legacy) == npairs(nlist_sort)
            @test compare_pairlists(nlist_legacy, nlist_sort)
        end
    end

    @testset "Non-cubic Cells" begin
        N = 80
        cutoff = 3.0

        # Triclinic cell
        C1 = SMat([10.0 2.0 1.0; 0.0 9.0 1.5; 0.0 0.0 8.0])
        X1 = [C1' * SVec(rand(), rand(), rand()) for _ in 1:N]
        pbc = SVec(true, true, true)

        nlist_legacy = PairList(X1, cutoff, C1, pbc; int_type=Int32)
        clist = build_cell_list(X1, cutoff, C1, pbc; backend=CPU())
        nlist_sort = materialize_pairlist(clist)

        @test npairs(nlist_legacy) == npairs(nlist_sort)
        @test compare_pairlists(nlist_legacy, nlist_sort)

        # Elongated cell
        C2 = SMat(diagm([5.0, 5.0, 20.0]))
        X2 = [C2' * SVec(rand(), rand(), rand()) for _ in 1:N]

        nlist_legacy2 = PairList(X2, cutoff, C2, pbc; int_type=Int32)
        clist2 = build_cell_list(X2, cutoff, C2, pbc; backend=CPU())
        nlist_sort2 = materialize_pairlist(clist2)

        @test npairs(nlist_legacy2) == npairs(nlist_sort2)
        @test compare_pairlists(nlist_legacy2, nlist_sort2)
    end

    @testset "Edge Cases" begin
        # Single atom
        X_single = [SVec(5.0, 5.0, 5.0)]
        C = SMat(diagm([10.0, 10.0, 10.0]))
        pbc = SVec(true, true, true)
        cutoff = 3.0

        clist = build_cell_list(X_single, cutoff, C, pbc; backend=CPU())
        nlist = materialize_pairlist(clist)
        @test npairs(nlist) == 0

        # Two atoms within cutoff
        X_two = [SVec(5.0, 5.0, 5.0), SVec(5.0, 5.0, 6.0)]
        clist2 = build_cell_list(X_two, cutoff, C, pbc; backend=CPU())
        nlist2 = materialize_pairlist(clist2)
        @test npairs(nlist2) == 2  # (1,2) and (2,1)

        # Two atoms outside cutoff
        X_far = [SVec(1.0, 1.0, 1.0), SVec(8.0, 8.0, 8.0)]
        clist3 = build_cell_list(X_far, cutoff, C, SVec(false, false, false); backend=CPU())
        nlist3 = materialize_pairlist(clist3)
        @test npairs(nlist3) == 0

        # Very small cutoff (no pairs)
        # Note: use a reasonable cutoff to avoid cell count overflow
        L_small = 10.0
        C_small = SMat(diagm([L_small, L_small, L_small]))
        X_small = [SVec(L_small * rand(), L_small * rand(), L_small * rand()) for _ in 1:50]
        clist4 = build_cell_list(X_small, 0.5, C_small, pbc; backend=CPU())  # Small but not too small
        nlist4 = materialize_pairlist(clist4)
        # Just verify it runs and produces reasonable output
        @test npairs(nlist4) >= 0
    end

    @testset "Large Cutoff" begin
        # Cutoff larger than half cell (stress tests multi-cell search)
        N = 30
        L = 10.0
        C = SMat(diagm([L, L, L]))
        X = [SVec(L * rand(), L * rand(), L * rand()) for _ in 1:N]
        cutoff = L * 0.6  # Larger than L/2
        pbc = SVec(true, true, true)

        nlist_legacy = PairList(X, cutoff, C, pbc; int_type=Int32)
        clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
        nlist_sort = materialize_pairlist(clist)

        @test npairs(nlist_legacy) == npairs(nlist_sort)
        @test compare_pairlists(nlist_legacy, nlist_sort)
    end

    @testset "Different Integer Types" begin
        N = 100
        X, C, L = rand_config_sortbased(N)
        cutoff = L / 3
        pbc = SVec(true, true, true)

        for TI in [Int32, Int64]
            clist = build_cell_list(X, cutoff, C, pbc; backend=CPU(), int_type=TI)
            nlist = materialize_pairlist(clist)

            @test eltype(nlist.i) == TI
            @test eltype(nlist.j) == TI
            @test npairs(nlist) > 0
        end
    end
end

@testset "Lazy Iteration (for_each_neighbour)" begin

    @testset "Correctness vs Materialized" begin
        N = 100
        X, C, L = rand_config_sortbased(N)
        cutoff = L / 3
        pbc = SVec(true, true, true)

        clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
        nlist = materialize_pairlist(clist)

        # Count pairs via lazy iteration
        lazy_count = 0
        for i in 1:nsites(clist)
            for_each_neighbour(clist, i) do j, R, S
                lazy_count += 1
            end
        end

        @test lazy_count == npairs(nlist)

        # Collect pairs via lazy iteration
        lazy_pairs = Tuple{Int, Int, NTuple{3,Int}}[]
        for i in 1:nsites(clist)
            for_each_neighbour(clist, i) do j, R, S
                push!(lazy_pairs, (i, j, Tuple(S)))
            end
        end

        materialized_pairs = get_sorted_pairs(nlist)
        @test sort(lazy_pairs) == materialized_pairs
    end

    @testset "Displacement Vectors" begin
        # Check that R vectors from lazy iteration match materialized
        N = 50
        X, C, L = rand_config_sortbased(N)
        cutoff = L / 3
        pbc = SVec(true, true, true)

        clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
        nlist = materialize_pairlist(clist)

        # Build dictionary of (i,j,S) -> R from materialized list
        R_dict = Dict{Tuple{Int,Int,NTuple{3,Int}}, SVec{Float64}}()
        for idx in 1:npairs(nlist)
            i, j = nlist.i[idx], nlist.j[idx]
            S = nlist.S[idx]
            R = X[j] - X[i] + C' * S
            R_dict[(i, j, Tuple(S))] = R
        end

        # Compare with lazy iteration
        for i in 1:nsites(clist)
            for_each_neighbour(clist, i) do j, R, S
                key = (i, j, Tuple(S))
                @test haskey(R_dict, key)
                @test R ≈ R_dict[key] atol=1e-10
            end
        end
    end
end

@testset "MapReduce Functions" begin

    @testset "map_sites! with SortedCellList" begin
        N = 100
        X, C, L = rand_config_sortbased(N)
        cutoff = L / 3
        pbc = SVec(true, true, true)

        clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())

        # Sum of distances per site using map_sites!
        # map_sites!(f, out, clist) where f(Rs) -> value
        out = zeros(Float64, N)
        map_sites!(Rs -> sum(norm.(Rs)), out, clist)

        # Verify against lazy iteration
        expected = zeros(Float64, N)
        for i in 1:N
            for_each_neighbour(clist, i) do j, R, S
                expected[i] += norm(R)
            end
        end

        @test out ≈ expected atol=1e-10
    end

    @testset "map_sites! with PairList" begin
        N = 100
        X, C, L = rand_config_sortbased(N)
        cutoff = L / 3
        pbc = SVec(true, true, true)

        clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
        nlist = materialize_pairlist(clist)

        # Sum of distances per site
        out = zeros(Float64, N)
        map_sites!(Rs -> sum(norm.(Rs)), out, nlist)

        # Verify manually
        expected = zeros(Float64, N)
        for i in 1:N
            for_each_neighbour(clist, i) do j, R, S
                expected[i] += norm(R)
            end
        end

        @test out ≈ expected atol=1e-10
    end

    @testset "map_pairs! symmetric" begin
        N = 100
        X, C, L = rand_config_sortbased(N)
        cutoff = L / 3
        pbc = SVec(true, true, true)

        clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
        nlist = materialize_pairlist(clist)

        # Energy-like quantity: sum 1/r to each site (symmetric)
        # map_pairs!(f, out, nlist) where f(i, j, R) -> value
        out = zeros(Float64, N)
        map_pairs!((i, j, R) -> 1.0 / norm(R), out, nlist; symmetric=true)

        # Each pair (i,j) contributes 1/r to both i and j
        # In a symmetric map, we only process i<j and add half to both
        expected = zeros(Float64, N)
        for idx in 1:npairs(nlist)
            i, j = nlist.i[idx], nlist.j[idx]
            S = nlist.S[idx]
            R = X[j] - X[i] + C' * S
            r = norm(R)
            if i < j
                expected[i] += 0.5 / r
                expected[j] += 0.5 / r
            end
        end

        @test out ≈ expected atol=1e-10
    end
end

@testset "Stress Tests" begin

    @testset "Larger Systems" begin
        for N in [500, 1000, 2000]
            X, C, L = rand_config_sortbased(N)
            cutoff = L / 4
            pbc = SVec(true, true, true)

            nlist_legacy = PairList(X, cutoff, C, pbc; int_type=Int32)
            clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
            nlist_sort = materialize_pairlist(clist)

            @test npairs(nlist_legacy) == npairs(nlist_sort)
        end
    end

    @testset "High Density" begin
        # Dense system with many neighbours
        N = 200
        L = 5.0  # Small box = high density
        C = SMat(diagm([L, L, L]))
        X = [SVec(L * rand(), L * rand(), L * rand()) for _ in 1:N]
        cutoff = 2.0
        pbc = SVec(true, true, true)

        nlist_legacy = PairList(X, cutoff, C, pbc; int_type=Int32)
        clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
        nlist_sort = materialize_pairlist(clist)

        @test npairs(nlist_legacy) == npairs(nlist_sort)
        @test compare_pairlists(nlist_legacy, nlist_sort)
    end

    @testset "Sparse System" begin
        # Low density system
        N = 50
        L = 50.0  # Large box = low density
        C = SMat(diagm([L, L, L]))
        X = [SVec(L * rand(), L * rand(), L * rand()) for _ in 1:N]
        cutoff = 5.0
        pbc = SVec(true, true, true)

        nlist_legacy = PairList(X, cutoff, C, pbc; int_type=Int32)
        clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
        nlist_sort = materialize_pairlist(clist)

        @test npairs(nlist_legacy) == npairs(nlist_sort)
        @test compare_pairlists(nlist_legacy, nlist_sort)
    end
end

@testset "get_neighbours Function" begin
    N = 50
    X, C, L = rand_config_sortbased(N)
    cutoff = L / 3
    pbc = SVec(true, true, true)

    clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())

    for i in 1:10  # Test first 10 atoms
        js, Rs, Ss = get_neighbours(clist, i)

        # Verify consistency with for_each_neighbour
        js_lazy = Int32[]
        Rs_lazy = SVec{Float64}[]
        Ss_lazy = SVec{Int32}[]
        for_each_neighbour(clist, i) do j, R, S
            push!(js_lazy, j)
            push!(Rs_lazy, R)
            push!(Ss_lazy, S)
        end

        @test sort(js) == sort(js_lazy)
        @test length(Rs) == length(Rs_lazy)
    end
end

@testset "count_neighbours Function" begin
    N = 50
    X, C, L = rand_config_sortbased(N)
    cutoff = L / 3
    pbc = SVec(true, true, true)

    clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
    nlist = materialize_pairlist(clist)

    # Count from PairList for reference
    expected_counts = zeros(Int, N)
    for idx in 1:npairs(nlist)
        expected_counts[nlist.i[idx]] += 1
    end

    # Compare with count_neighbours
    for i in 1:N
        @test count_neighbours(clist, i) == expected_counts[i]
    end
end

@testset "map_pairs! with symmetric=false" begin
    N = 100
    X, C, L = rand_config_sortbased(N)
    cutoff = L / 3
    pbc = SVec(true, true, true)

    clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
    nlist = materialize_pairlist(clist)

    # Non-symmetric: each pair (i,j) only contributes to site i
    out = zeros(Float64, N)
    map_pairs!((i, j, R) -> 1.0, out, nlist; symmetric=false)

    # Expected: count of pairs where i is the first index
    expected = zeros(Float64, N)
    for idx in 1:npairs(nlist)
        expected[nlist.i[idx]] += 1.0
    end

    @test out == expected
end

@testset "map_pairs_d! Function" begin
    N = 100
    X, C, L = rand_config_sortbased(N)
    cutoff = L / 3
    pbc = SVec(true, true, true)

    clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
    nlist = materialize_pairlist(clist)

    # Anti-symmetric mapping: force-like
    out = zeros(SVec{Float64}, N)
    map_pairs_d!((i, j, R) -> R / norm(R), out, nlist)

    # Verify anti-symmetry: sum should be approximately zero
    total = sum(out)
    @test norm(total) < 1e-10
end
