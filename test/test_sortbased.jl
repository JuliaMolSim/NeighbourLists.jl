# Tests for sort-based cell list implementation
# Uses shared utilities from test_utils.jl (included by runtests.jl)

using NeighbourLists: SortedCellList, map_sites!, map_pairs!, map_pairs_d!,
                      count_neighbours, neighbours, CPU

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

@testset "neighbours / count_neighbours cell list" begin
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
        js, Rs, Ss = neighbours(clist, i)
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

# Regression tests for issue #6 — atoms outside the unit cell with PBC.
#
# `_build_sorted_celllist` bins atoms by their coordinates wrapped to
# `[0, L)`, but stores the *unwrapped* positions in `clist.X` /
# `clist.X_orig`. `for_each_neighbour` then computes
# `R = X[j] - X[i] + cell' * cell_shift` from unwrapped positions, so the
# pair vector misses the minimum image whenever the source atom started
# outside `[0, L)`. Symptom: empty neighbour lists for atoms placed
# outside the box, even when their periodic-image distance is well below
# cutoff.

# Collect (sorted-key, distance) pairs from a SortedCellList by lazy
# iteration. Uses an order-canonicalised key so the result is
# permutation-invariant in atom labelling.
function _collect_pair_multiset(clist)
    pairs = Tuple{Int, Int, Float64}[]
    for i in 1:nsites(clist)
        for_each_neighbour(clist, i) do j, R, S
            push!(pairs, (min(i, Int(j)), max(i, Int(j)),
                          round(sqrt(sum(abs2, R)); digits = 10)))
        end
    end
    return sort(pairs)
end

@testset "for_each_neighbour with atoms outside [0, L) (#6)" begin

    @testset "two atoms 1 Å apart along y" begin
        L = 8.0
        C = cubic_cell(L)
        cutoff = 1.5
        pbc = FULL_PBC

        # Reference: both atoms inside [0, L)
        X_in = [SVec(0.5, 0.5, 0.5), SVec(0.5, 1.5, 0.5)]
        # Bug-triggering configuration: atom 2 just outside [0, L)
        # along y. Distance |X_in[1] - X_in[2]| == |X_bug[1] - X_bug[2]|
        # under the periodic image, so the pair geometry is identical.
        X_bug = [SVec(0.5, 0.5, 0.5), SVec(0.5, -0.5, 0.5)]
        # Atom 2 displaced by a full period (-L) from its X_in counterpart
        # (1.5 - L = -6.5). Wraps onto X_in[2].
        X_far = [SVec(0.5, 0.5, 0.5), SVec(0.5, -6.5, 0.5)]

        clist_in  = build_cell_list(X_in,  cutoff, C, pbc; backend = CPU())
        clist_bug = build_cell_list(X_bug, cutoff, C, pbc; backend = CPU())
        clist_far = build_cell_list(X_far, cutoff, C, pbc; backend = CPU())

        ref = _collect_pair_multiset(clist_in)
        @test length(ref) == 2                       # ordered pairs (1,2),(2,1)
        @test _collect_pair_multiset(clist_bug) == ref
        @test _collect_pair_multiset(clist_far) == ref
    end

    @testset "individually shifting some atoms by ±L on a periodic axis is a no-op" begin
        # Uniform translation doesn't trigger the bug (every atom shifts
        # cells consistently). The bug shows when *some* atoms are
        # inside `[0, L)` and others are outside — exactly the situation
        # H2O builders create. We perturb each atom by a deterministic
        # integer multiple of `L` along each periodic axis (open axes
        # are left alone), cycling through {-1, 0, +1} so the resulting
        # configuration mixes in-box and out-of-box atoms.
        for pbc in all_pbc_cases()
            X, C, L = rand_config(60)
            cutoff = L * 0.25
            ref = _collect_pair_multiset(
                build_cell_list(X, cutoff, C, pbc; backend = CPU()))

            X_perturbed = map(enumerate(X)) do (i, x)
                # Deterministic per-axis shift in {-1, 0, +1}, distinct
                # per (axis, atom) so different atoms land in different
                # periodic images.
                shifts = (mod(i,     3) - 1,
                          mod(i + 1, 3) - 1,
                          mod(i + 2, 3) - 1)
                offset = SVec(ntuple(
                    α -> pbc[α] ? L * Float64(shifts[α]) : 0.0, 3)...)
                x + offset
            end
            got = _collect_pair_multiset(
                build_cell_list(X_perturbed, cutoff, C, pbc; backend = CPU()))
            @test got == ref
        end
    end
end
