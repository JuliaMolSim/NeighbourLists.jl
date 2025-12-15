#!/usr/bin/env julia
"""
Performance benchmark for NeighbourLists.jl

Compares legacy (linked-list), sort-based CPU (multi-threaded), and GPU implementations.
Results are used in the README.md benchmark table.

Usage:
    # Multi-threaded CPU benchmark (16 threads) + GPU
    julia --project -t 16 scripts/benchmark.jl

    # Single-threaded CPU benchmark
    julia --project -t 1 scripts/benchmark.jl
"""

using Pkg
Pkg.activate(dirname(@__DIR__))

using NeighbourLists
using NeighbourLists: build_cell_list, materialize_pairlist, npairs, PairList, 
                      SVec, CPU 
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

# Check for PrettyTables
try
    using PrettyTables
catch
    println("Installing PrettyTables...")
    Pkg.add("PrettyTables")
    using PrettyTables
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

    # Collect benchmark results
    results = []

    for n_atoms in n_atoms_list
        X, cell, pbc, L = generate_system(n_atoms, density)

        # Get pair count
        clist = build_cell_list(X, cutoff, cell, pbc; backend=CPU())
        nlist = materialize_pairlist(clist)
        n_pairs = npairs(nlist)

        row = Dict(
            :atoms => n_atoms,
            :pairs => format_pairs(n_pairs),
            :legacy_time => "-",
            :cpu_time => "-",
            :gpu_time => "-",
            :speedup_vs_legacy => "-"
        )

        # Convert to SVec for legacy API
        X_svec = [SVec(x...) for x in X]

        # Legacy (linked-list) benchmark - suppress deprecation warnings
        if n_atoms <= 100000
            t_legacy = @belapsed begin
                nlist = PairList($X_svec, $cutoff, $cell, $pbc; int_type=Int32)
            end samples=5 evals=1
            row[:legacy_time] = format_time(t_legacy)
            row[:t_legacy_raw] = t_legacy
        end

        # Sort-based CPU benchmark
        if n_atoms <= 100000
            t_cpu = @belapsed begin
                clist = build_cell_list($X, $cutoff, $cell, $pbc; backend=CPU())
                nlist = materialize_pairlist(clist)
            end samples=5 evals=1
            row[:cpu_time] = format_time(t_cpu)
            row[:t_cpu_raw] = t_cpu
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

            row[:gpu_time] = format_time(t_gpu)
            row[:t_gpu_raw] = t_gpu

            # Speedup vs legacy
            if haskey(row, :t_legacy_raw)
                speedup = row[:t_legacy_raw] / t_gpu
                row[:speedup_vs_legacy] = @sprintf("%.1fx", speedup)
            end
        end

        push!(results, row)
    end

    # Build table data
    atoms_col = [r[:atoms] for r in results]
    pairs_col = [r[:pairs] for r in results]
    legacy_col = [r[:legacy_time] for r in results]
    cpu_col = [r[:cpu_time] for r in results]

    if cuda_available
        gpu_col = [r[:gpu_time] for r in results]
        speedup_col = [r[:speedup_vs_legacy] for r in results]

        data = hcat(atoms_col, pairs_col, legacy_col, cpu_col, gpu_col, speedup_col)
        header = ["Atoms", "Pairs", "Legacy", "CPU ($(Threads.nthreads())T)", "GPU", "Speedup"]
    else
        data = hcat(atoms_col, pairs_col, legacy_col, cpu_col)
        header = ["Atoms", "Pairs", "Legacy", "CPU ($(Threads.nthreads())T)"]
    end

    # Print table using PrettyTables (v3.x API)
    pretty_table(data;
        backend = :markdown,
        column_labels = header,
        alignment = [:r for _ in 1:length(header)]
    )

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
