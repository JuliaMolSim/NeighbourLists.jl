# Tests validating NeighbourLists.jl against Python matscipy neighbourlist
# Uses PythonCall to interface with matscipy

using NeighbourLists
using NeighbourLists: SMat, SVec, nsites, npairs, build_cell_list, materialize_pairlist
using Test
using LinearAlgebra
using StaticArrays

# Check for PythonCall and CondaPkg availability
pythoncall_available = false
try
    using CondaPkg
    CondaPkg.resolve()  # Ensure Python environment is set up
    using PythonCall
    global pythoncall_available = true
catch e
    @info "PythonCall/CondaPkg not available, skipping matscipy tests: $e"
end

if pythoncall_available
    # Import matscipy
    matscipy_available = false
    local matscipy_neighbours_mod = nothing
    local numpy_mod = nothing
    try
        matscipy_neighbours_mod = pyimport("matscipy.neighbours")
        numpy_mod = pyimport("numpy")
        global matscipy_available = true
    catch e
        @info "matscipy not available in Python environment: $e"
    end

    if matscipy_available
        @info "matscipy available, running comparison tests"

        # Store references for use in function
        const MATSCIPY_NEIGHBOURS = matscipy_neighbours_mod
        const NUMPY = numpy_mod

        """
        Get neighbour list from matscipy for comparison.
        Returns (i, j, S) arrays where S is the shift vector.
        """
        function matscipy_neighbourlist(positions, cell, pbc, cutoff)
            # Convert Julia arrays to numpy arrays
            # Flatten positions to 1D array, then reshape
            pos_flat = collect(reinterpret(Float64, positions))
            pos_np = NUMPY.array(pos_flat).reshape(length(positions), 3)

            # Cell matrix - transpose for row-major format expected by Python
            cell_flat = collect(transpose(cell))[:]
            cell_np = NUMPY.array(cell_flat).reshape(3, 3)

            # PBC array
            pbc_np = NUMPY.array(collect(pbc))

            # Call matscipy neighbour_list
            # Returns i, j, S (shift vectors)
            result = MATSCIPY_NEIGHBOURS.neighbour_list("ijS", positions=pos_np,
                                                         cell=cell_np, pbc=pbc_np,
                                                         cutoff=Float64(cutoff))

            i_py, j_py, S_py = result

            # Convert back to Julia arrays (0-indexed to 1-indexed)
            i_jl = pyconvert(Vector{Int}, i_py) .+ 1
            j_jl = pyconvert(Vector{Int}, j_py) .+ 1
            S_np = pyconvert(Matrix{Int}, S_py)
            S_jl = [SVec{Int}(S_np[k, 1], S_np[k, 2], S_np[k, 3]) for k in 1:size(S_np, 1)]

            return i_jl, j_jl, S_jl
        end

        """
        Compare Julia neighbour list with matscipy result.
        Returns true if they match (order-independent comparison).
        """
        function compare_with_matscipy(nlist, i_ms, j_ms, S_ms)
            # Build set of (i, j, S) tuples from Julia
            jl_pairs = Set{Tuple{Int, Int, Tuple{Int,Int,Int}}}()
            for idx in 1:npairs(nlist)
                push!(jl_pairs, (nlist.i[idx], nlist.j[idx], Tuple(nlist.S[idx])))
            end

            # Build set from matscipy
            ms_pairs = Set{Tuple{Int, Int, Tuple{Int,Int,Int}}}()
            for idx in 1:length(i_ms)
                push!(ms_pairs, (i_ms[idx], j_ms[idx], Tuple(S_ms[idx])))
            end

            return jl_pairs == ms_pairs
        end

        # Helper to generate random configurations
        function rand_config_matscipy(N; density=0.05)
            volume = N / density
            L = volume^(1/3)
            C = SMat(diagm([L, L, L]))
            X = [SVec(L * rand(), L * rand(), L * rand()) for _ in 1:N]
            return X, C, L
        end

        @testset "Matscipy Comparison" begin

            @testset "Basic Cubic Cell - Full PBC" begin
                N = 100
                X, C, L = rand_config_matscipy(N)
                cutoff = L / 3
                pbc = SVec(true, true, true)

                # Julia implementation
                clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
                nlist = materialize_pairlist(clist)

                # Matscipy implementation
                i_ms, j_ms, S_ms = matscipy_neighbourlist(X, C, pbc, cutoff)

                @test npairs(nlist) == length(i_ms)
                @test compare_with_matscipy(nlist, i_ms, j_ms, S_ms)
            end

            @testset "No PBC" begin
                N = 80
                X, C, L = rand_config_matscipy(N)
                cutoff = L / 4
                pbc = SVec(false, false, false)

                clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
                nlist = materialize_pairlist(clist)

                i_ms, j_ms, S_ms = matscipy_neighbourlist(X, C, pbc, cutoff)

                @test npairs(nlist) == length(i_ms)
                @test compare_with_matscipy(nlist, i_ms, j_ms, S_ms)
            end

            @testset "Mixed PBC" begin
                N = 80
                X, C, L = rand_config_matscipy(N)
                cutoff = L / 4

                # Test several mixed PBC combinations
                pbc_cases = [
                    SVec(true, false, false),
                    SVec(false, true, false),
                    SVec(false, false, true),
                    SVec(true, true, false),
                    SVec(true, false, true),
                    SVec(false, true, true),
                ]

                for pbc in pbc_cases
                    clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
                    nlist = materialize_pairlist(clist)

                    i_ms, j_ms, S_ms = matscipy_neighbourlist(X, C, pbc, cutoff)

                    @test npairs(nlist) == length(i_ms)
                    @test compare_with_matscipy(nlist, i_ms, j_ms, S_ms)
                end
            end

            @testset "Non-cubic Cells" begin
                N = 60
                cutoff = 3.0
                pbc = SVec(true, true, true)

                # Triclinic cell
                C1 = SMat([10.0 2.0 1.0; 0.0 9.0 1.5; 0.0 0.0 8.0])
                X1 = [C1' * SVec(rand(), rand(), rand()) for _ in 1:N]

                clist1 = build_cell_list(X1, cutoff, C1, pbc; backend=CPU())
                nlist1 = materialize_pairlist(clist1)

                i_ms1, j_ms1, S_ms1 = matscipy_neighbourlist(X1, C1, pbc, cutoff)

                @test npairs(nlist1) == length(i_ms1)
                @test compare_with_matscipy(nlist1, i_ms1, j_ms1, S_ms1)

                # Elongated cell
                C2 = SMat(diagm([5.0, 5.0, 20.0]))
                X2 = [C2' * SVec(rand(), rand(), rand()) for _ in 1:N]

                clist2 = build_cell_list(X2, cutoff, C2, pbc; backend=CPU())
                nlist2 = materialize_pairlist(clist2)

                i_ms2, j_ms2, S_ms2 = matscipy_neighbourlist(X2, C2, pbc, cutoff)

                @test npairs(nlist2) == length(i_ms2)
                @test compare_with_matscipy(nlist2, i_ms2, j_ms2, S_ms2)
            end

            @testset "Multiple Random Configurations" begin
                for trial in 1:10
                    N = rand(30:150)
                    X, C, L = rand_config_matscipy(N)
                    cutoff = L / 4 + rand() * L / 4
                    pbc = SVec(rand(Bool), rand(Bool), rand(Bool))

                    clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
                    nlist = materialize_pairlist(clist)

                    i_ms, j_ms, S_ms = matscipy_neighbourlist(X, C, pbc, cutoff)

                    @test npairs(nlist) == length(i_ms)
                    @test compare_with_matscipy(nlist, i_ms, j_ms, S_ms)
                end
            end

            @testset "Edge Cases" begin
                pbc = SVec(true, true, true)

                # Single atom (no pairs)
                X_single = [SVec(5.0, 5.0, 5.0)]
                C = SMat(diagm([10.0, 10.0, 10.0]))
                cutoff = 3.0

                clist = build_cell_list(X_single, cutoff, C, pbc; backend=CPU())
                nlist = materialize_pairlist(clist)
                i_ms, j_ms, S_ms = matscipy_neighbourlist(X_single, C, pbc, cutoff)

                @test npairs(nlist) == length(i_ms) == 0

                # Two atoms within cutoff
                X_two = [SVec(5.0, 5.0, 5.0), SVec(5.0, 5.0, 6.0)]

                clist2 = build_cell_list(X_two, cutoff, C, pbc; backend=CPU())
                nlist2 = materialize_pairlist(clist2)
                i_ms2, j_ms2, S_ms2 = matscipy_neighbourlist(X_two, C, pbc, cutoff)

                @test npairs(nlist2) == length(i_ms2)
                @test compare_with_matscipy(nlist2, i_ms2, j_ms2, S_ms2)
            end

            @testset "Larger Systems" begin
                for N in [200, 500, 1000]
                    X, C, L = rand_config_matscipy(N)
                    cutoff = L / 4
                    pbc = SVec(true, true, true)

                    clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
                    nlist = materialize_pairlist(clist)

                    i_ms, j_ms, S_ms = matscipy_neighbourlist(X, C, pbc, cutoff)

                    @test npairs(nlist) == length(i_ms)
                    # For larger systems, just check counts match to save time
                    if N <= 500
                        @test compare_with_matscipy(nlist, i_ms, j_ms, S_ms)
                    end
                end
            end

            @testset "Legacy PairList vs Matscipy" begin
                # Also validate the legacy implementation
                N = 100
                X, C, L = rand_config_matscipy(N)
                cutoff = L / 3
                pbc = SVec(true, true, true)

                # Legacy implementation
                nlist_legacy = PairList(X, cutoff, C, pbc; int_type=Int32)

                # Matscipy
                i_ms, j_ms, S_ms = matscipy_neighbourlist(X, C, pbc, cutoff)

                @test npairs(nlist_legacy) == length(i_ms)
                @test compare_with_matscipy(nlist_legacy, i_ms, j_ms, S_ms)
            end

        end

    else
        @testset "Matscipy Comparison (Skipped - matscipy not available)" begin
            @test_skip "matscipy not available"
        end
    end

else
    @testset "Matscipy Comparison (Skipped - PythonCall not available)" begin
        @test_skip "PythonCall not available"
    end
end
