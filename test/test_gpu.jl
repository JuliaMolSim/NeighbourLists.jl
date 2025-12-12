# GPU tests for NeighbourLists.jl
# Uses shared utilities from test_utils.jl (included by runtests.jl)

# Use cuda_available from runtests.jl (set via check_cuda_available())
if cuda_available
    using CUDA
    @info "CUDA available: $(CUDA.name(CUDA.device()))"

    @testset "GPU Cell List" begin

        @testset "Basic Construction" begin
            X, C, L = rand_config(100)
            X_gpu = CuArray(X)
            clist = build_cell_list(X_gpu, L/3, C, SVec(true,true,true))

            @test clist isa SortedCellList
            @test nsites(clist) == 100
            @test clist.X_orig isa CuArray
            @test clist.perm isa CuArray
        end

        @testset "Pair Materialization" begin
            X, C, L = rand_config(100)
            X_gpu = CuArray(X)
            nlist = materialize_pairlist(build_cell_list(X_gpu, L/3, C, SVec(true,true,true)))

            @test nlist.i isa CuArray
            @test nlist.j isa CuArray
            @test nlist.S isa CuArray
            @test npairs(nlist) > 0
        end
    end

    @testset "GPU vs CPU Correctness" begin

        @testset "Random Configs" begin
            for _ in 1:5
                test_cpu_vs_gpu(CuArray; N=rand(50:200))
            end
        end

        @testset "All PBC Combinations" begin
            test_all_pbc_cpu_vs_gpu(CuArray)
        end

        @testset "Triclinic Cell" begin
            test_triclinic_cell() do X, C, cutoff, pbc
                nlist_cpu = materialize_pairlist(build_cell_list(X, cutoff, C, pbc; backend=CPU()))
                nlist_gpu = materialize_pairlist(build_cell_list(CuArray(X), cutoff, C, pbc))
                @test npairs(nlist_cpu) == npairs(nlist_gpu)
                @test compare_cpu_gpu_full(nlist_cpu, nlist_gpu)
            end
        end

        @testset "Full (i,j,S) Match" begin
            X, C, L = rand_config(100)
            nlist_cpu = materialize_pairlist(build_cell_list(X, L/3, C, SVec(true,true,true); backend=CPU()))
            nlist_gpu = materialize_pairlist(build_cell_list(CuArray(X), L/3, C, SVec(true,true,true)))
            @test compare_cpu_gpu_full(nlist_cpu, nlist_gpu)
        end
    end

    @testset "GPU Edge Cases" begin
        C = cubic_cell(10.0)
        pbc = SVec(true, true, true)

        # Single atom
        nlist = materialize_pairlist(build_cell_list(CuArray([SVec(5.0,5.0,5.0)]), 3.0, C, pbc))
        @test npairs(nlist) == 0

        # Two atoms
        X = [SVec(5.0, 5.0, 5.0), SVec(5.0, 5.0, 6.0)]
        nlist_cpu = materialize_pairlist(build_cell_list(X, 3.0, C, pbc; backend=CPU()))
        nlist_gpu = materialize_pairlist(build_cell_list(CuArray(X), 3.0, C, pbc))
        @test npairs(nlist_cpu) == npairs(nlist_gpu)
    end

    @testset "GPU Large Systems" begin
        test_sizes(sizes=[500, 1000, 2000]) do X, C, L, cutoff
            nlist_cpu = materialize_pairlist(build_cell_list(X, cutoff, C, SVec(true,true,true); backend=CPU()))
            nlist_gpu = materialize_pairlist(build_cell_list(CuArray(X), cutoff, C, SVec(true,true,true)))
            @test npairs(nlist_cpu) == npairs(nlist_gpu)
        end
    end

    @testset "GPU High Density" begin
        L = 5.0
        C = cubic_cell(L)
        X = [SVec(L * rand(), L * rand(), L * rand()) for _ in 1:200]
        nlist_cpu = materialize_pairlist(build_cell_list(X, 2.0, C, SVec(true,true,true); backend=CPU()))
        nlist_gpu = materialize_pairlist(build_cell_list(CuArray(X), 2.0, C, SVec(true,true,true)))
        @test npairs(nlist_cpu) == npairs(nlist_gpu)
        @test compare_cpu_gpu_full(nlist_cpu, nlist_gpu)
    end

else
    @testset "GPU Tests (Skipped)" begin
        @test_skip "CUDA not available"
    end
end
