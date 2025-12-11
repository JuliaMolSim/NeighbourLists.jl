# GPU tests for NeighbourLists.jl
# These tests are only run if CUDA is available and functional

using NeighbourLists
using NeighbourLists: SMat, SVec, SortedCellList, nsites, npairs,
                      build_cell_list, materialize_pairlist
using Test
using LinearAlgebra
using StaticArrays

# Check for CUDA availability
cuda_available = false
try
    using CUDA
    global cuda_available = CUDA.functional()
catch
    @info "CUDA not available, skipping GPU tests"
end

if cuda_available
    @info "CUDA available: $(CUDA.name(CUDA.device()))"

    # Helper to generate random configurations
    function rand_config_gpu(N; density=0.05)
        volume = N / density
        L = volume^(1/3)
        C = SMat(diagm([L, L, L]))
        X = [SVec(L * rand(), L * rand(), L * rand()) for _ in 1:N]
        return X, C, L
    end

    # Helper to compare CPU and GPU pair lists
    function compare_cpu_gpu_pairs(nlist_cpu, nlist_gpu)
        cpu_pairs = sort(collect(zip(nlist_cpu.i, nlist_cpu.j)))
        gpu_pairs = sort(collect(zip(Array(nlist_gpu.i), Array(nlist_gpu.j))))
        return cpu_pairs == gpu_pairs
    end

    # Helper to compare shift vectors
    function compare_cpu_gpu_shifts(nlist_cpu, nlist_gpu)
        cpu_shifts = sort([Tuple(s) for s in nlist_cpu.S])
        gpu_shifts = sort([Tuple(s) for s in Array(nlist_gpu.S)])
        return cpu_shifts == gpu_shifts
    end

    @testset "GPU Cell List Construction" begin

        @testset "Basic GPU Construction" begin
            N = 100
            X, C, L = rand_config_gpu(N)
            cutoff = L / 3
            pbc = SVec(true, true, true)

            # CPU version
            clist_cpu = build_cell_list(X, cutoff, C, pbc; backend=CPU())

            # GPU version
            X_gpu = CuArray(X)
            clist_gpu = build_cell_list(X_gpu, cutoff, C, pbc)

            @test clist_gpu isa SortedCellList
            @test nsites(clist_gpu) == N
            @test clist_gpu.cutoff == cutoff

            # Check that arrays are on GPU
            @test clist_gpu.X_orig isa CuArray
            @test clist_gpu.perm isa CuArray
            @test clist_gpu.cell_offsets isa CuArray
        end

        @testset "GPU Pair Materialization" begin
            N = 100
            X, C, L = rand_config_gpu(N)
            cutoff = L / 3
            pbc = SVec(true, true, true)

            # CPU version
            clist_cpu = build_cell_list(X, cutoff, C, pbc; backend=CPU())
            nlist_cpu = materialize_pairlist(clist_cpu)

            # GPU version
            X_gpu = CuArray(X)
            clist_gpu = build_cell_list(X_gpu, cutoff, C, pbc)
            nlist_gpu = materialize_pairlist(clist_gpu)

            # Check pair arrays are on GPU
            @test nlist_gpu.i isa CuArray
            @test nlist_gpu.j isa CuArray
            @test nlist_gpu.S isa CuArray

            # Compare pair counts
            @test npairs(nlist_cpu) == npairs(nlist_gpu)

            # Compare actual pairs
            @test compare_cpu_gpu_pairs(nlist_cpu, nlist_gpu)
        end
    end

    @testset "GPU vs CPU Correctness" begin

        @testset "Multiple Random Configurations" begin
            for trial in 1:5
                N = rand(50:200)
                X, C, L = rand_config_gpu(N)
                cutoff = L / 4 + rand() * L / 4
                pbc = SVec(rand(Bool), rand(Bool), rand(Bool))

                # CPU
                clist_cpu = build_cell_list(X, cutoff, C, pbc; backend=CPU())
                nlist_cpu = materialize_pairlist(clist_cpu)

                # GPU
                X_gpu = CuArray(X)
                clist_gpu = build_cell_list(X_gpu, cutoff, C, pbc)
                nlist_gpu = materialize_pairlist(clist_gpu)

                @test npairs(nlist_cpu) == npairs(nlist_gpu)
                @test compare_cpu_gpu_pairs(nlist_cpu, nlist_gpu)
            end
        end

        @testset "All PBC Combinations" begin
            N = 80
            X, C, L = rand_config_gpu(N)
            cutoff = L / 3

            pbc_cases = [
                SVec(true, true, true),
                SVec(false, false, false),
                SVec(true, false, false),
                SVec(false, true, false),
                SVec(false, false, true),
                SVec(true, true, false),
                SVec(true, false, true),
                SVec(false, true, true),
            ]

            for pbc in pbc_cases
                clist_cpu = build_cell_list(X, cutoff, C, pbc; backend=CPU())
                nlist_cpu = materialize_pairlist(clist_cpu)

                X_gpu = CuArray(X)
                clist_gpu = build_cell_list(X_gpu, cutoff, C, pbc)
                nlist_gpu = materialize_pairlist(clist_gpu)

                @test npairs(nlist_cpu) == npairs(nlist_gpu)
                @test compare_cpu_gpu_pairs(nlist_cpu, nlist_gpu)
            end
        end

        @testset "Non-cubic Cells" begin
            N = 60
            cutoff = 3.0
            pbc = SVec(true, true, true)

            # Triclinic
            C1 = SMat([10.0 2.0 1.0; 0.0 9.0 1.5; 0.0 0.0 8.0])
            X1 = [C1' * SVec(rand(), rand(), rand()) for _ in 1:N]

            clist_cpu = build_cell_list(X1, cutoff, C1, pbc; backend=CPU())
            nlist_cpu = materialize_pairlist(clist_cpu)

            X1_gpu = CuArray(X1)
            clist_gpu = build_cell_list(X1_gpu, cutoff, C1, pbc)
            nlist_gpu = materialize_pairlist(clist_gpu)

            @test npairs(nlist_cpu) == npairs(nlist_gpu)
            @test compare_cpu_gpu_pairs(nlist_cpu, nlist_gpu)
        end

        @testset "Shift Vectors Match" begin
            N = 100
            X, C, L = rand_config_gpu(N)
            cutoff = L / 3
            pbc = SVec(true, true, true)

            clist_cpu = build_cell_list(X, cutoff, C, pbc; backend=CPU())
            nlist_cpu = materialize_pairlist(clist_cpu)

            X_gpu = CuArray(X)
            clist_gpu = build_cell_list(X_gpu, cutoff, C, pbc)
            nlist_gpu = materialize_pairlist(clist_gpu)

            # Build full (i,j,S) tuples and compare
            cpu_full = sort(collect(zip(nlist_cpu.i, nlist_cpu.j, [Tuple(s) for s in nlist_cpu.S])))
            gpu_full = sort(collect(zip(Array(nlist_gpu.i), Array(nlist_gpu.j), [Tuple(s) for s in Array(nlist_gpu.S)])))

            @test cpu_full == gpu_full
        end
    end

    @testset "GPU Edge Cases" begin

        @testset "Single Atom" begin
            X = [SVec(5.0, 5.0, 5.0)]
            C = SMat(diagm([10.0, 10.0, 10.0]))
            pbc = SVec(true, true, true)
            cutoff = 3.0

            X_gpu = CuArray(X)
            clist_gpu = build_cell_list(X_gpu, cutoff, C, pbc)
            nlist_gpu = materialize_pairlist(clist_gpu)

            @test npairs(nlist_gpu) == 0
        end

        @testset "Two Atoms" begin
            X = [SVec(5.0, 5.0, 5.0), SVec(5.0, 5.0, 6.0)]
            C = SMat(diagm([10.0, 10.0, 10.0]))
            pbc = SVec(true, true, true)
            cutoff = 3.0

            clist_cpu = build_cell_list(X, cutoff, C, pbc; backend=CPU())
            nlist_cpu = materialize_pairlist(clist_cpu)

            X_gpu = CuArray(X)
            clist_gpu = build_cell_list(X_gpu, cutoff, C, pbc)
            nlist_gpu = materialize_pairlist(clist_gpu)

            @test npairs(nlist_cpu) == npairs(nlist_gpu)
        end

        @testset "No Pairs (small cutoff)" begin
            N = 50
            X, C, L = rand_config_gpu(N)
            cutoff = 0.01  # Very small
            pbc = SVec(true, true, true)

            X_gpu = CuArray(X)
            clist_gpu = build_cell_list(X_gpu, cutoff, C, pbc)
            nlist_gpu = materialize_pairlist(clist_gpu)

            @test npairs(nlist_gpu) == 0
        end
    end

    @testset "GPU Larger Systems" begin
        for N in [500, 1000, 2000]
            X, C, L = rand_config_gpu(N)
            cutoff = L / 4
            pbc = SVec(true, true, true)

            clist_cpu = build_cell_list(X, cutoff, C, pbc; backend=CPU())
            nlist_cpu = materialize_pairlist(clist_cpu)

            X_gpu = CuArray(X)
            clist_gpu = build_cell_list(X_gpu, cutoff, C, pbc)
            nlist_gpu = materialize_pairlist(clist_gpu)

            @test npairs(nlist_cpu) == npairs(nlist_gpu)

            # For larger systems, just verify counts match (full comparison is slow)
            if N <= 1000
                @test compare_cpu_gpu_pairs(nlist_cpu, nlist_gpu)
            end
        end
    end

    @testset "GPU High Density" begin
        N = 200
        L = 5.0  # Small box
        C = SMat(diagm([L, L, L]))
        X = [SVec(L * rand(), L * rand(), L * rand()) for _ in 1:N]
        cutoff = 2.0
        pbc = SVec(true, true, true)

        clist_cpu = build_cell_list(X, cutoff, C, pbc; backend=CPU())
        nlist_cpu = materialize_pairlist(clist_cpu)

        X_gpu = CuArray(X)
        clist_gpu = build_cell_list(X_gpu, cutoff, C, pbc)
        nlist_gpu = materialize_pairlist(clist_gpu)

        @test npairs(nlist_cpu) == npairs(nlist_gpu)
        @test compare_cpu_gpu_pairs(nlist_cpu, nlist_gpu)
    end

else
    @testset "GPU Tests (Skipped - CUDA not available)" begin
        @test_skip "CUDA not available"
    end
end
