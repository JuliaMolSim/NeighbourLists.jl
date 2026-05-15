# GPU tests for NeighbourLists.jl
# Uses shared utilities from test_utils.jl (included by runtests.jl)
# Supports multiple GPU backends: CUDA, AMDGPU, Metal

using NeighbourLists: CPU
using KernelAbstractions
using KernelAbstractions: @kernel, @index

if gpu_available()
    backend_name = gpu_backend()
    @info "Running GPU tests with $backend_name backend"

    @testset "GPU Cell List ($backend_name)" begin

        @testset "Basic Construction" begin
            X, C, L = rand_config(100)
            X_gpu = to_gpu_array(X)
            clist = build_cell_list(X_gpu, L/3, C, FULL_PBC)

            @test clist isa SortedCellList
            @test nsites(clist) == 100
            @test clist.X_orig isa gpu_array_type()
            @test clist.perm isa gpu_array_type()
        end

        @testset "Pair Materialization" begin
            X, C, L = rand_config(100)
            X_gpu = to_gpu_array(X)
            nlist = materialize_pairlist(build_cell_list(X_gpu, L/3, C, FULL_PBC))

            @test nlist.i isa gpu_array_type()
            @test nlist.j isa gpu_array_type()
            @test nlist.S isa gpu_array_type()
            @test npairs(nlist) > 0
        end
    end

    @testset "GPU vs CPU Correctness ($backend_name)" begin

        @testset "Random Configs" begin
            for _ in 1:5
                test_cpu_vs_gpu(to_gpu_array; N=rand(50:200))
            end
        end

        @testset "All PBC Combinations" begin
            test_all_pbc_cpu_vs_gpu(to_gpu_array)
        end

        @testset "Triclinic Cell" begin
            test_triclinic_cell() do X, C, cutoff, pbc
                nlist_cpu = materialize_pairlist(build_cell_list(X, cutoff, C, pbc; backend=CPU()))
                nlist_gpu = materialize_pairlist(build_cell_list(to_gpu_array(X), cutoff, C, pbc))
                @test npairs(nlist_cpu) == npairs(nlist_gpu)
                @test compare_cpu_gpu_full(nlist_cpu, nlist_gpu)
            end
        end

        @testset "Full (i,j,S) Match" begin
            X, C, L = rand_config(100)
            nlist_cpu = materialize_pairlist(build_cell_list(X, L/3, C, FULL_PBC; backend=CPU()))
            nlist_gpu = materialize_pairlist(build_cell_list(to_gpu_array(X), L/3, C, FULL_PBC))
            @test compare_cpu_gpu_full(nlist_cpu, nlist_gpu)
        end
    end

    @testset "GPU Edge Cases ($backend_name)" begin
        test_edge_cases_cpu_vs_gpu(to_gpu_array)
    end

    @testset "GPU Large Systems ($backend_name)" begin
        test_large_systems_cpu_vs_gpu(to_gpu_array)
    end

    @testset "GPU High Density ($backend_name)" begin
        L = 5.0
        C = cubic_cell(L)
        X = [SVec(L * rand(), L * rand(), L * rand()) for _ in 1:200]
        nlist_cpu = materialize_pairlist(build_cell_list(X, 2.0, C, FULL_PBC; backend=CPU()))
        nlist_gpu = materialize_pairlist(build_cell_list(to_gpu_array(X), 2.0, C, FULL_PBC))
        @test npairs(nlist_cpu) == npairs(nlist_gpu)
        @test compare_cpu_gpu_full(nlist_cpu, nlist_gpu)
    end

    # Regression for #37: passing a SortedCellList as a kernel argument
    # fails GPU compilation when no `Adapt.adapt_structure` rule is
    # defined on it — the inner CuArray fields don't get recursively
    # adapted to `CuDeviceArray`. Internal NL.jl kernels work around
    # this by unpacking the fields manually; external consumers calling
    # `for_each_neighbour` inside a custom `@kernel` cannot.
    @testset "SortedCellList in user kernels ($backend_name)" begin
        X, C, L = rand_config(50)
        X_gpu = to_gpu_array(X)
        clist = build_cell_list(X_gpu, L/3, C, FULL_PBC)

        @kernel function _count_neighbours_via_clist!(out, c)
            i = @index(Global, Linear)
            n = Ref(zero(eltype(out)))
            for_each_neighbour(c, i) do j, R, S
                n[] += one(eltype(out))
            end
            out[i] = n[]
        end

        backend = KernelAbstractions.get_backend(X_gpu)
        out = KernelAbstractions.zeros(backend, Int32, nsites(clist))
        _count_neighbours_via_clist!(backend)(out, clist; ndrange = nsites(clist))
        KernelAbstractions.synchronize(backend)

        # Compare against a materialised-pairlist neighbour-count reference.
        nlist = materialize_pairlist(clist)
        expected = zeros(Int32, nsites(clist))
        for i in Array(nlist.i)
            expected[i] += one(Int32)
        end
        @test Array(out) == expected
    end

else
    @testset "GPU Tests (Skipped)" begin
        @test_skip "No GPU backend available (CUDA, AMDGPU, or Metal)"
    end
end
