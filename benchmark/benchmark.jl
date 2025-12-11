#!/usr/bin/env julia
"""
Performance benchmark for NeighbourLists.jl

Compares:
1. Legacy serial code (linked-list based)
2. New sort-based CPU (single-threaded)
3. New sort-based CPU (multi-threaded)
4. GPU (CUDA) if available
"""

using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(path=dirname(@__DIR__))
Pkg.instantiate()

using NeighbourLists
using StaticArrays
using BenchmarkTools
using Printf
using LinearAlgebra

# Check for CUDA
cuda_available = false
try
    using CUDA
    global cuda_available = CUDA.functional()
    if cuda_available
        println("CUDA available: ", CUDA.name(CUDA.device()))
    end
catch
    println("CUDA not available")
end

# Generate random atomic positions in a cubic cell
function generate_system(n_atoms::Int, density::Float64=0.05)
    # density in atoms per Å³
    volume = n_atoms / density
    L = volume^(1/3)

    X = [SVector{3,Float64}(L * rand(), L * rand(), L * rand()) for _ in 1:n_atoms]
    cell = SMatrix{3,3,Float64,9}(L, 0, 0, 0, L, 0, 0, 0, L)
    pbc = SVector{3,Bool}(true, true, true)

    return X, cell, pbc, L
end

# Benchmark functions
function bench_legacy(X, cutoff, cell, pbc; samples=10)
    # Legacy linked-list based construction
    times = Float64[]
    for _ in 1:samples
        t = @elapsed begin
            nlist = PairList(Vector(X), cutoff, Matrix(cell), Tuple(pbc);
                            int_type=Int32, fixcell=false)
        end
        push!(times, t)
    end
    return times
end

function bench_sortbased_cpu(X, cutoff, cell, pbc; samples=10)
    # New sort-based construction (CPU)
    times = Float64[]
    for _ in 1:samples
        t = @elapsed begin
            clist = build_cell_list(X, cutoff, cell, pbc; backend=CPU())
            nlist = materialize_pairlist(clist)
        end
        push!(times, t)
    end
    return times
end

function bench_sortbased_cpu_celllist_only(X, cutoff, cell, pbc; samples=10)
    # Just cell list construction (no pair materialization)
    times = Float64[]
    for _ in 1:samples
        t = @elapsed begin
            clist = build_cell_list(X, cutoff, cell, pbc; backend=CPU())
        end
        push!(times, t)
    end
    return times
end

function bench_lazy_iteration(X, cutoff, cell, pbc; samples=10)
    # Build cell list once, then time lazy iteration
    clist = build_cell_list(X, cutoff, cell, pbc; backend=CPU())

    times = Float64[]
    for _ in 1:samples
        t = @elapsed begin
            total_pairs = 0
            for i in 1:nsites(clist)
                for_each_neighbour(clist, i) do j, R, S
                    total_pairs += 1
                end
            end
        end
        push!(times, t)
    end
    return times
end

function format_time(t_seconds)
    if t_seconds < 1e-3
        return @sprintf("%.2f μs", t_seconds * 1e6)
    elseif t_seconds < 1
        return @sprintf("%.2f ms", t_seconds * 1e3)
    else
        return @sprintf("%.2f s", t_seconds)
    end
end

function run_benchmark(n_atoms_list, cutoff, density)
    println("\n" * "="^70)
    println("NeighbourLists.jl Performance Benchmark")
    println("="^70)
    println("Cutoff: $(cutoff) Å, Density: $(density) atoms/Å³")
    println("="^70)

    results = Dict{Int, Dict{String, Float64}}()

    for n_atoms in n_atoms_list
        println("\n--- $n_atoms atoms ---")
        X, cell, pbc, L = generate_system(n_atoms, density)
        println("Box size: $(round(L, digits=2)) Å")

        # Warmup
        _ = PairList(Vector(X), cutoff, Matrix(cell), Tuple(pbc); int_type=Int32, fixcell=false)
        _ = build_cell_list(X, cutoff, cell, pbc; backend=CPU())

        results[n_atoms] = Dict{String, Float64}()

        # Legacy serial
        print("  Legacy (linked-list): ")
        times = bench_legacy(X, cutoff, cell, pbc; samples=5)
        t_legacy = minimum(times)
        results[n_atoms]["legacy"] = t_legacy
        println(format_time(t_legacy))

        # Sort-based CPU (cell list only)
        print("  Sort-based cell list: ")
        times = bench_sortbased_cpu_celllist_only(X, cutoff, cell, pbc; samples=5)
        t_celllist = minimum(times)
        results[n_atoms]["celllist"] = t_celllist
        println(format_time(t_celllist))

        # Sort-based CPU (full)
        print("  Sort-based full:      ")
        times = bench_sortbased_cpu(X, cutoff, cell, pbc; samples=5)
        t_sortbased = minimum(times)
        results[n_atoms]["sortbased"] = t_sortbased
        println(format_time(t_sortbased))

        # Lazy iteration
        print("  Lazy iteration:       ")
        times = bench_lazy_iteration(X, cutoff, cell, pbc; samples=5)
        t_lazy = minimum(times)
        results[n_atoms]["lazy"] = t_lazy
        println(format_time(t_lazy))

        # Count pairs for reference
        nlist = PairList(Vector(X), cutoff, Matrix(cell), Tuple(pbc); int_type=Int32, fixcell=false)
        println("  Total pairs: $(npairs(nlist))")

        # Speedup
        speedup = t_legacy / t_sortbased
        println("  Speedup (sort vs legacy): $(round(speedup, digits=2))x")
    end

    return results
end

function run_gpu_benchmark(n_atoms, cutoff, density)
    if !cuda_available
        println("\nGPU benchmark skipped (CUDA not available)")
        return
    end

    println("\n" * "="^70)
    println("GPU Benchmark ($(n_atoms) atoms)")
    println("="^70)

    X, cell, pbc, L = generate_system(n_atoms, density)
    println("Box size: $(round(L, digits=2)) Å, Cutoff: $(cutoff) Å")

    # CPU baseline
    print("  CPU cell list build: ")
    t_cpu = @belapsed build_cell_list($X, $cutoff, $cell, $pbc; backend=CPU()) samples=5
    println(format_time(t_cpu))

    # CPU materialization
    clist_cpu = build_cell_list(X, cutoff, cell, pbc; backend=CPU())
    print("  CPU materialize:     ")
    t_mat_cpu = @belapsed materialize_pairlist($clist_cpu) samples=5
    println(format_time(t_mat_cpu))

    nlist = materialize_pairlist(clist_cpu)
    println("  Total pairs: $(npairs(nlist))")

    # Note: Full GPU implementation requires GPU-side pair materialization
    # which would need more kernel development
    println("\n  Note: Full GPU pair computation requires additional kernel development")
    println("  The sort-based architecture is GPU-ready for future optimization")
end

function run_threading_benchmark(n_atoms, cutoff, density)
    println("\n" * "="^70)
    println("Threading Benchmark ($(n_atoms) atoms)")
    println("="^70)

    X, cell, pbc, L = generate_system(n_atoms, density)
    println("Box size: $(round(L, digits=2)) Å, Cutoff: $(cutoff) Å")
    println("Available threads: $(Threads.nthreads())")

    # Warmup
    clist = build_cell_list(X, cutoff, cell, pbc; backend=CPU())
    nlist = materialize_pairlist(clist)
    println("Total pairs: $(npairs(nlist))")

    # Test map_sites! with threading (uses @threads automatically)
    println("\nmap_sites! benchmark (computing sum of distances per atom):")
    out = zeros(Float64, length(X))

    print("  With @threads: ")
    t_threaded = @belapsed begin
        fill!($out, 0.0)
        map_sites!(Rs -> sum(norm.(Rs)), $out, $nlist)
    end samples=5
    println(format_time(t_threaded))

    # Test lazy iteration with threading
    println("\nLazy parallel iteration (counting pairs per atom):")
    counts = zeros(Int, length(X))

    print("  Serial:    ")
    t_serial = @belapsed begin
        fill!($counts, 0)
        for i in 1:length($X)
            for_each_neighbour($clist, i) do j, R, S
                $counts[i] += 1
            end
        end
    end samples=5
    println(format_time(t_serial))

    print("  Threaded:  ")
    t_parallel = @belapsed begin
        fill!($counts, 0)
        Threads.@threads for i in 1:length($X)
            for_each_neighbour($clist, i) do j, R, S
                $counts[i] += 1
            end
        end
    end samples=5
    speedup = t_serial / t_parallel
    println(format_time(t_parallel), " ($(round(speedup, digits=2))x speedup)")
end

# Main
println("Julia threads: ", Threads.nthreads())

# Small to medium system sizes
n_atoms_list = [100, 500, 1000, 5000, 10000]
cutoff = 5.0  # Å
density = 0.05  # atoms/Å³ (roughly liquid density)

results = run_benchmark(n_atoms_list, cutoff, density)

# Threading benchmark for larger system
if Threads.nthreads() > 1
    run_threading_benchmark(10000, cutoff, density)
else
    println("\nNote: Run with multiple threads (julia -t N) for threading benchmark")
end

# GPU benchmark
run_gpu_benchmark(50000, cutoff, density)

println("\n" * "="^70)
println("Benchmark complete")
println("="^70)
