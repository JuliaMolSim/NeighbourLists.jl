#!/usr/bin/env julia
"""
Performance benchmark for NeighbourLists.jl

Compares sort-based implementation on CPU (single-threaded and multi-threaded) vs GPU.
Results are used in the README.md benchmark table.

Usage:
    # Multi-threaded CPU benchmark (16 threads) + GPU
    julia --project -t 16 scripts/benchmark.jl

    # Single-threaded CPU benchmark (for 1T column in README)
    julia --project -t 1 scripts/benchmark.jl

To get complete benchmark data for README:
1. Run with -t 1 to get CPU (1T) times
2. Run with -t 16 (or desired thread count) to get CPU (nT) and GPU times
"""

using Pkg
Pkg.activate(dirname(@__DIR__))

using NeighbourLists
using NeighbourLists: build_cell_list, materialize_pairlist, npairs
using StaticArrays
using LinearAlgebra
using Printf

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

# Check for BenchmarkTools
try
    using BenchmarkTools
catch
    println("Installing BenchmarkTools...")
    Pkg.add("BenchmarkTools")
    using BenchmarkTools
end

# Generate random atomic positions in a cubic cell
function generate_system(n_atoms::Int, density::Float64=0.05)
    volume = n_atoms / density
    L = volume^(1/3)
    X = [SVector{3,Float64}(L * rand(), L * rand(), L * rand()) for _ in 1:n_atoms]
    cell = SMatrix{3,3,Float64,9}(L, 0, 0, 0, L, 0, 0, 0, L)
    pbc = SVector{3,Bool}(true, true, true)
    return X, cell, pbc, L
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

function format_pairs(n)
    if n < 1000
        return string(n)
    elseif n < 1_000_000
        return @sprintf("%.0fk", n / 1000)
    else
        return @sprintf("%.1fM", n / 1_000_000)
    end
end

function run_benchmark(; n_atoms_list=[1000, 5000, 10000, 50000, 100000],
                        cutoff=5.0, density=0.05)
    println("\n" * "="^70)
    println("NeighbourLists.jl Performance Benchmark")
    println("="^70)
    println("Julia version: ", VERSION)
    println("Julia threads: ", Threads.nthreads())
    if cuda_available
        println("GPU: ", CUDA.name(CUDA.device()))
    end
    println("Cutoff: $(cutoff) Å, Density: $(density) atoms/Å³")
    println()

    # Build header
    header = "| Atoms | Pairs | CPU ($(Threads.nthreads())T) |"
    separator = "|------:|------:|----------:|"
    if cuda_available
        header *= "      GPU | Speedup |"
        separator *= "---------:|--------:|"
    end

    println(header)
    println(separator)

    for n_atoms in n_atoms_list
        X, cell, pbc, L = generate_system(n_atoms, density)

        # Get pair count
        clist = build_cell_list(X, cutoff, cell, pbc; backend=CPU())
        nlist = materialize_pairlist(clist)
        n_pairs = npairs(nlist)

        row = "| $(lpad(n_atoms, 6)) | $(lpad(format_pairs(n_pairs), 5)) |"

        # CPU benchmark
        if n_atoms <= 100000
            t_cpu = @belapsed begin
                clist = build_cell_list($X, $cutoff, $cell, $pbc; backend=CPU())
                nlist = materialize_pairlist(clist)
            end samples=5 evals=1
            row *= " $(lpad(format_time(t_cpu), 9)) |"
        else
            t_cpu = nothing
            row *= " $(lpad("-", 9)) |"
        end

        # GPU benchmark
        if cuda_available
            X_gpu = CuArray(X)
            # Warmup
            clist_gpu = build_cell_list(X_gpu, cutoff, cell, pbc)
            nlist_gpu = materialize_pairlist(clist_gpu)

            t_gpu = @belapsed begin
                clist = build_cell_list($X_gpu, $cutoff, $cell, $pbc)
                nlist = materialize_pairlist(clist)
            end samples=5 evals=1

            row *= " $(lpad(format_time(t_gpu), 8)) |"

            # Speedup
            if t_cpu !== nothing
                speedup = t_cpu / t_gpu
                row *= " $(lpad(@sprintf("%.1fx", speedup), 7)) |"
            else
                row *= " $(lpad("-", 7)) |"
            end
        end

        println(row)
    end

    # GPU throughput for large system
    if cuda_available
        println()
        max_atoms = maximum(n_atoms_list)
        X, cell, pbc, L = generate_system(max_atoms, density)
        X_gpu = CuArray(X)
        clist_gpu = build_cell_list(X_gpu, cutoff, cell, pbc)
        nlist_gpu = materialize_pairlist(clist_gpu)

        t_gpu = @belapsed begin
            clist = build_cell_list($X_gpu, $cutoff, $cell, $pbc)
            nlist = materialize_pairlist(clist)
        end samples=5 evals=1

        throughput = npairs(nlist_gpu) / t_gpu / 1e6
        println("GPU throughput: $(round(Int, throughput)) million pairs/second")
    end
end

# Main
println("NeighbourLists.jl Benchmark")

cutoff = 5.0  # Å
density = 0.05  # atoms/Å³ (roughly liquid density)

# Run benchmarks
run_benchmark(; cutoff=cutoff, density=density)

println("\n" * "="^70)
println("Benchmark complete")
println("="^70)
