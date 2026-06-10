# GPU tests for NeighbourLists.jl
# Uses shared utilities from test_utils.jl (included by runtests.jl)
# Supports multiple GPU backends: CUDA, AMDGPU, Metal

using NeighbourLists: CPU
using KernelAbstractions
using KernelAbstractions: @kernel, @index

if gpu_available()
    backend_name = gpu_backend()
    eltypes = gpu_supported_eltypes()
    @info "Running GPU tests with $backend_name backend ($(join(eltypes, ", ")))"

    for T in eltypes

        @testset "GPU Cell List ($backend_name, $T)" begin

            @testset "Basic Construction" begin
                X, C, L = rand_config(100; T=T)
                X_gpu = to_gpu_array(X)
                clist = build_cell_list(X_gpu, L/3, C, FULL_PBC)

                @test clist isa SortedCellList
                @test nsites(clist) == 100
                @test clist.X_orig isa gpu_array_type()
                @test clist.perm isa gpu_array_type()
            end

            @testset "Pair Materialization" begin
                X, C, L = rand_config(100; T=T)
                X_gpu = to_gpu_array(X)
                nlist = materialize_pairlist(build_cell_list(X_gpu, L/3, C, FULL_PBC))

                @test nlist.i isa gpu_array_type()
                @test nlist.j isa gpu_array_type()
                @test nlist.S isa gpu_array_type()
                @test npairs(nlist) > 0
            end
        end

        @testset "GPU vs CPU Correctness ($backend_name, $T)" begin

            @testset "Random Configs" begin
                for _ in 1:5
                    test_cpu_vs_gpu(to_gpu_array; N=rand(50:200), T=T)
                end
            end

            @testset "All PBC Combinations" begin
                test_all_pbc_cpu_vs_gpu(to_gpu_array; T=T)
            end

            @testset "Triclinic Cell" begin
                test_triclinic_cell(; T=T) do X, C, cutoff, pbc
                    nlist_cpu = materialize_pairlist(build_cell_list(X, cutoff, C, pbc; backend=CPU()))
                    nlist_gpu = materialize_pairlist(build_cell_list(to_gpu_array(X), cutoff, C, pbc))
                    @test npairs(nlist_cpu) == npairs(nlist_gpu)
                    @test compare_cpu_gpu_full(nlist_cpu, nlist_gpu)
                end
            end

            @testset "Full (i,j,S) Match" begin
                X, C, L = rand_config(100; T=T)
                nlist_cpu = materialize_pairlist(build_cell_list(X, L/3, C, FULL_PBC; backend=CPU()))
                nlist_gpu = materialize_pairlist(build_cell_list(to_gpu_array(X), L/3, C, FULL_PBC))
                @test compare_cpu_gpu_full(nlist_cpu, nlist_gpu)
            end
        end

        @testset "GPU Edge Cases ($backend_name, $T)" begin
            test_edge_cases_cpu_vs_gpu(to_gpu_array; T=T)
        end

        @testset "GPU Large Systems ($backend_name, $T)" begin
            if backend_name == :Metal
                # AcceleratedKernels.sortperm! is broken on Metal for
                # n ≳ 512 (silently wrong results, possible GPU hang) —
                # see issue #42. Cap the system size until that is fixed.
                test_large_systems_cpu_vs_gpu(to_gpu_array; sizes=[500], T=T)
                @test_skip "sizes > 512 disabled on Metal (#42)"
            else
                test_large_systems_cpu_vs_gpu(to_gpu_array; T=T)
            end
        end

        @testset "GPU High Density ($backend_name, $T)" begin
            L = T(5)
            C = cubic_cell(L)
            X = [SVec{T}(L * rand(T), L * rand(T), L * rand(T)) for _ in 1:200]
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
        @testset "SortedCellList in user kernels ($backend_name, $T)" begin
            X, C, L = rand_config(50; T=T)
            X_gpu = to_gpu_array(X)
            clist = build_cell_list(X_gpu, L/3, C, FULL_PBC)

            # NB: accumulate directly into `out` rather than a local `Ref`:
            # `Ref` heap-allocates inside the kernel, which CUDA tolerates
            # but Metal cannot compile.
            @kernel function _count_neighbours_via_clist!(out, c)
                i = @index(Global, Linear)
                out[i] = zero(eltype(out))
                for_each_neighbour(c, i) do j, R, S
                    out[i] += one(eltype(out))
                end
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

    end # for T in eltypes

    if gpu_backend() == :Metal
        @testset "GPU Float64 (Skipped on Metal)" begin
            @test_skip "Metal does not support Float64"
        end
    end

else
    @testset "GPU Tests (Skipped)" begin
        @test_skip "No GPU backend available (CUDA, AMDGPU, or Metal)"
    end
end
