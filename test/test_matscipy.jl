# Tests validating NeighbourLists.jl against Python matscipy neighbourlist
# Uses shared utilities from test_utils.jl (included by runtests.jl)

pythoncall_available = check_pythoncall_available()

if !pythoncall_available
    @testset "Matscipy Comparison (Skipped)" begin
        @test_skip "PythonCall not available"
    end
else
    using CondaPkg, PythonCall

    matscipy_available = false
    local MATSCIPY, NUMPY
    try
        MATSCIPY = pyimport("matscipy.neighbours")
        NUMPY = pyimport("numpy")
        global matscipy_available = true
    catch e
        @info "matscipy not available: $e"
    end

    if !matscipy_available
        @testset "Matscipy Comparison (Skipped)" begin
            @test_skip "matscipy not available"
        end
    else
        @info "matscipy available, running comparison tests"

        # Matscipy interface
        function matscipy_neighbourlist(X, C, pbc, cutoff)
            pos_np = NUMPY.array(collect(reinterpret(Float64, X))).reshape(length(X), 3)
            cell_np = NUMPY.array(collect(transpose(C))[:]).reshape(3, 3)
            pbc_np = NUMPY.array(collect(pbc))

            i_py, j_py, S_py = MATSCIPY.neighbour_list("ijS", positions=pos_np,
                                                        cell=cell_np, pbc=pbc_np,
                                                        cutoff=Float64(cutoff))

            i_jl = pyconvert(Vector{Int}, i_py) .+ 1
            j_jl = pyconvert(Vector{Int}, j_py) .+ 1
            S_np = pyconvert(Matrix{Int}, S_py)
            S_jl = [SVec{Int}(S_np[k, 1], S_np[k, 2], S_np[k, 3]) for k in 1:size(S_np, 1)]
            return i_jl, j_jl, S_jl
        end

        function compare_with_matscipy(nlist, i_ms, j_ms, S_ms)
            jl_set = pairs_to_set(nlist.i, nlist.j, nlist.S)
            ms_set = pairs_to_set(i_ms, j_ms, S_ms)
            return jl_set == ms_set
        end

        # Parametric test helper
        function test_vs_matscipy(X, C, cutoff, pbc)
            nlist = materialize_pairlist(build_cell_list(X, cutoff, C, pbc; backend=CPU()))
            i_ms, j_ms, S_ms = matscipy_neighbourlist(X, C, pbc, cutoff)
            @test npairs(nlist) == length(i_ms)
            @test compare_with_matscipy(nlist, i_ms, j_ms, S_ms)
        end

        @testset "Matscipy Comparison" begin

            @testset "Basic Cubic - Full PBC" begin
                X, C, L = rand_config(100)
                test_vs_matscipy(X, C, L/3, FULL_PBC)
            end

            @testset "No PBC" begin
                X, C, L = rand_config(80)
                test_vs_matscipy(X, C, L/4, NO_PBC)
            end

            @testset "All PBC Combinations" begin
                X, C, L = rand_config(80)
                for pbc in all_pbc_cases()
                    test_vs_matscipy(X, C, L/4, pbc)
                end
            end

            @testset "Triclinic Cell" begin
                test_triclinic_cell() do X, C, cutoff, pbc
                    test_vs_matscipy(X, C, cutoff, pbc)
                end
            end

            @testset "Random Configurations" begin
                for _ in 1:10
                    X, C, L = rand_config(rand(30:150))
                    cutoff = L * (0.25 + 0.25 * rand())
                    pbc = SVec(rand(Bool), rand(Bool), rand(Bool))
                    test_vs_matscipy(X, C, cutoff, pbc)
                end
            end

            @testset "Edge Cases" begin
                C = cubic_cell(10.0)

                # Single atom
                test_vs_matscipy([SVec(5.0, 5.0, 5.0)], C, 3.0, FULL_PBC)

                # Two atoms
                test_vs_matscipy([SVec(5.0, 5.0, 5.0), SVec(5.0, 5.0, 6.0)], C, 3.0, FULL_PBC)
            end

            @testset "Large Systems" begin
                for N in [200, 500, 1000]
                    X, C, L = rand_config(N)
                    nlist = materialize_pairlist(build_cell_list(X, L/4, C, FULL_PBC; backend=CPU()))
                    i_ms, _, _ = matscipy_neighbourlist(X, C, FULL_PBC, L/4)
                    @test npairs(nlist) == length(i_ms)
                end
            end

            @testset "Legacy PairList vs Matscipy" begin
                X, C, L = rand_config(100)
                nlist = PairList(X, L/3, C, FULL_PBC; int_type=Int32)
                i_ms, j_ms, S_ms = matscipy_neighbourlist(X, C, FULL_PBC, L/3)
                @test npairs(nlist) == length(i_ms)
                @test compare_with_matscipy(nlist, i_ms, j_ms, S_ms)
            end

            if gpu_available()
                backend_name = gpu_backend()
                @info "Running GPU ($backend_name) vs matscipy validation"

                @testset "GPU vs Matscipy ($backend_name)" begin
                    X, C, L = rand_config(100)
                    nlist = materialize_pairlist(build_cell_list(to_gpu_array(X), L/3, C, FULL_PBC))
                    i_ms, j_ms, S_ms = matscipy_neighbourlist(X, C, FULL_PBC, L/3)

                    @test npairs(nlist) == length(i_ms)
                    gpu_set = pairs_to_set(Array(nlist.i), Array(nlist.j), Array(nlist.S))
                    ms_set = pairs_to_set(i_ms, j_ms, S_ms)
                    @test gpu_set == ms_set
                end

                @testset "GPU Multiple Configs ($backend_name)" begin
                    for _ in 1:5
                        X, C, L = rand_config(rand(50:150))
                        cutoff = L * (0.25 + 0.25 * rand())
                        pbc = SVec(rand(Bool), rand(Bool), rand(Bool))
                        nlist = materialize_pairlist(build_cell_list(to_gpu_array(X), cutoff, C, pbc))
                        i_ms, _, _ = matscipy_neighbourlist(X, C, pbc, cutoff)
                        @test npairs(nlist) == length(i_ms)
                    end
                end
            else
                @testset "GPU vs Matscipy (Skipped)" begin
                    @test_skip "No GPU backend available"
                end
            end
        end
    end
end
