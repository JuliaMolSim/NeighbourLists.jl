# Shared test utilities for NeighbourLists.jl
# This module consolidates common test helpers to reduce code duplication

using NeighbourLists
using NeighbourLists: SMat, SVec, nsites, npairs, SortedCellList,
                      build_cell_list, materialize_pairlist, for_each_neighbour,
                      neigs, count_neighbours, CPU 
using LinearAlgebra
using StaticArrays
using Test

# ==================== Constants ====================

"""Full periodic boundary conditions (all directions)"""
const FULL_PBC = SVec(true, true, true)

"""No periodic boundary conditions"""
const NO_PBC = SVec(false, false, false)

# ==================== Configuration Generators ====================

"""
    rand_config(N; density=0.05, T=Float64)

Generate a random cubic configuration with N atoms and element type `T`.
Returns (X, C, L) where X is positions, C is cell matrix, L is box length.
"""
function rand_config(N; density=0.05, T=Float64)
    volume = N / density
    L = T(volume^(1/3))
    C = SMat{T}(diagm([L, L, L]))
    X = [SVec{T}(L * rand(T), L * rand(T), L * rand(T)) for _ in 1:N]
    return X, C, L
end

"""
    all_pbc_cases()

Return all 8 combinations of periodic boundary conditions for 3D.
"""
function all_pbc_cases()
    return [
        SVec(true, true, true),    # Full PBC
        SVec(false, false, false), # No PBC
        SVec(true, false, false),  # PBC in x only
        SVec(false, true, false),  # PBC in y only
        SVec(false, false, true),  # PBC in z only
        SVec(true, true, false),   # PBC in x,y
        SVec(true, false, true),   # PBC in x,z
        SVec(false, true, true),   # PBC in y,z
    ]
end

"""
    cubic_cell(L)

Create a cubic cell matrix with side length L.
"""
cubic_cell(L) = SMat(diagm([L, L, L]))

# ==================== Pair List Comparison ====================

"""
    get_sorted_pairs(nlist)

Extract pairs from a PairList as sorted vector of (i, j, S) tuples.
"""
function get_sorted_pairs(nlist)
    return sort(collect(zip(nlist.i, nlist.j, [Tuple(s) for s in nlist.S])))
end

"""
    pairs_to_set(i, j, S)

Convert pair arrays (i, j, S) to a Set of tuples for comparison.
"""
function pairs_to_set(i, j, S)
    result = Set{Tuple{Int, Int, Tuple{Int,Int,Int}}}()
    for idx in eachindex(i)
        push!(result, (i[idx], j[idx], Tuple(S[idx])))
    end
    return result
end

"""
    compare_pairlists(nlist1, nlist2)

Compare two PairLists for equality (order-independent).
Returns true if they contain the same pairs.
"""
function compare_pairlists(nlist1, nlist2)
    pairs1 = get_sorted_pairs(nlist1)
    pairs2 = get_sorted_pairs(nlist2)
    return pairs1 == pairs2
end

# ==================== GPU/CPU Comparison ====================

"""
    compare_cpu_gpu_pairs(nlist_cpu, nlist_gpu)

Compare CPU and GPU pair lists for equality.
GPU arrays are copied to CPU for comparison.
"""
function compare_cpu_gpu_pairs(nlist_cpu, nlist_gpu)
    cpu_pairs = sort(collect(zip(nlist_cpu.i, nlist_cpu.j)))
    gpu_pairs = sort(collect(zip(Array(nlist_gpu.i), Array(nlist_gpu.j))))
    return cpu_pairs == gpu_pairs
end

"""
    compare_cpu_gpu_shifts(nlist_cpu, nlist_gpu)

Compare shift vectors between CPU and GPU pair lists.
"""
function compare_cpu_gpu_shifts(nlist_cpu, nlist_gpu)
    cpu_shifts = sort([Tuple(s) for s in nlist_cpu.S])
    gpu_shifts = sort([Tuple(s) for s in Array(nlist_gpu.S)])
    return cpu_shifts == gpu_shifts
end

"""
    compare_cpu_gpu_full(nlist_cpu, nlist_gpu)

Full comparison of CPU and GPU pair lists including shifts.
"""
function compare_cpu_gpu_full(nlist_cpu, nlist_gpu)
    cpu_set = pairs_to_set(nlist_cpu.i, nlist_cpu.j, nlist_cpu.S)
    gpu_set = pairs_to_set(Array(nlist_gpu.i), Array(nlist_gpu.j), Array(nlist_gpu.S))
    return cpu_set == gpu_set
end

# ==================== Backend Availability ====================

# GPU backend state (initialized once)
const GPU_STATE = Dict{Symbol, Any}(
    :initialized => false,
    :backend => :none,
    :to_gpu => identity,
    :gpu_array_type => nothing,
)

"""
    init_gpu_backend!()

Initialize GPU backend detection. Tries CUDA, then AMDGPU, then Metal.
Sets global GPU state that can be queried with `gpu_available()`, `gpu_backend()`, etc.
"""
function init_gpu_backend!()
    GPU_STATE[:initialized] && return GPU_STATE[:backend]

    # Try CUDA first
    try
        CUDA = Base.require(Main, :CUDA)
        # Use invokelatest to handle world age issues with dynamically loaded modules
        if Base.invokelatest(CUDA.functional)
            device_name = Base.invokelatest(CUDA.name, Base.invokelatest(CUDA.device))
            @info "GPU Backend: CUDA" device=device_name
            GPU_STATE[:backend] = :CUDA
            GPU_STATE[:to_gpu] = x -> Base.invokelatest(CUDA.cu, x)
            GPU_STATE[:gpu_array_type] = CUDA.CuArray
            GPU_STATE[:initialized] = true
            return :CUDA
        else
            @info "CUDA loaded but not functional"
        end
    catch e
        @info "CUDA not available" exception=e
    end

    # Try AMDGPU
    try
        AMDGPU = Base.require(Main, :AMDGPU)
        if Base.invokelatest(AMDGPU.functional)
            @info "GPU Backend: AMDGPU"
            GPU_STATE[:backend] = :AMDGPU
            GPU_STATE[:to_gpu] = x -> Base.invokelatest(AMDGPU.roc, x)
            GPU_STATE[:gpu_array_type] = AMDGPU.ROCArray
            GPU_STATE[:initialized] = true
            return :AMDGPU
        else
            @info "AMDGPU loaded but not functional"
        end
    catch e
        @info "AMDGPU not available" exception=e
    end

    # Try Metal
    try
        Metal = Base.require(Main, :Metal)
        if Base.invokelatest(Metal.functional)
            @info "GPU Backend: Metal"
            GPU_STATE[:backend] = :Metal
            GPU_STATE[:to_gpu] = x -> Base.invokelatest(Metal.mtl, x)
            GPU_STATE[:gpu_array_type] = Metal.MtlArray
            GPU_STATE[:initialized] = true
            return :Metal
        else
            @info "Metal loaded but not functional"
        end
    catch e
        @info "Metal not available" exception=e
    end

    @info "No GPU backend available, using CPU"
    GPU_STATE[:backend] = :none
    GPU_STATE[:to_gpu] = identity
    GPU_STATE[:gpu_array_type] = nothing
    GPU_STATE[:initialized] = true
    return :none
end

"""
    gpu_available()

Returns true if a GPU backend is available.
"""
function gpu_available()
    GPU_STATE[:initialized] || init_gpu_backend!()
    return GPU_STATE[:backend] != :none
end

"""
    gpu_backend()

Returns the name of the available GPU backend as a Symbol (:CUDA, :AMDGPU, :Metal, or :none).
"""
function gpu_backend()
    GPU_STATE[:initialized] || init_gpu_backend!()
    return GPU_STATE[:backend]
end

"""
    gpu_supported_eltypes()

Element types supported by the available GPU backend. Metal (Apple GPUs)
does not support Float64, so only Float32 is tested there.
"""
gpu_supported_eltypes() = gpu_backend() == :Metal ? (Float32,) : (Float64, Float32)

"""
    to_gpu_array(x)

Convert array `x` to the appropriate GPU array type for the available backend.
Returns `x` unchanged if no GPU is available.
"""
function to_gpu_array(x)
    GPU_STATE[:initialized] || init_gpu_backend!()
    return GPU_STATE[:to_gpu](x)
end

"""
    gpu_array_type()

Returns the GPU array type for the available backend (e.g., CuArray, ROCArray, MtlArray),
or `nothing` if no GPU is available.
"""
function gpu_array_type()
    GPU_STATE[:initialized] || init_gpu_backend!()
    return GPU_STATE[:gpu_array_type]
end

# Legacy function for backwards compatibility
"""
    check_cuda_available()

DEPRECATED: Use `gpu_available()` instead.
Check if CUDA is available and functional.
"""
function check_cuda_available()
    GPU_STATE[:initialized] || init_gpu_backend!()
    return GPU_STATE[:backend] == :CUDA
end

"""
    check_pythoncall_available()

Check if PythonCall and CondaPkg are available.
"""
function check_pythoncall_available()
    try
        Base.require(Main, :CondaPkg)
        Base.require(Main, :PythonCall)
        return true
    catch
        return false
    end
end

# ==================== Test Helpers ====================

"""
    rand_non_cubic_cell(T=Float64)

Generate a random non-cubic (triclinic) cell matrix for testing.
"""
function rand_non_cubic_cell(T=Float64)
    return SMat{T}([10.0 2.0 1.0; 0.0 9.0 1.5; 0.0 0.0 8.0])
end

"""
    rand_config_triclinic(N; cell=rand_non_cubic_cell())

Generate random configuration in a non-cubic cell.
"""
function rand_config_triclinic(N; cell=rand_non_cubic_cell())
    T = eltype(cell)
    X = [cell' * SVec{T}(rand(T), rand(T), rand(T)) for _ in 1:N]
    return X, cell
end

# ==================== High-Level Parametric Test Runners ====================

"""
    test_legacy_vs_sortbased(; N=100, density=0.05, cutoff_frac=0.25, pbc=SVec(true,true,true))

Test that sort-based implementation matches legacy linked-list implementation.
"""
function test_legacy_vs_sortbased(; N=100, density=0.05, cutoff_frac=0.25,
                                    pbc=SVec(true, true, true))
    X, C, L = rand_config(N; density=density)
    cutoff = L * cutoff_frac

    nlist_legacy = PairList(X, cutoff, C, pbc; int_type=Int32)
    clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
    nlist_sort = materialize_pairlist(clist)

    @test npairs(nlist_legacy) == npairs(nlist_sort)
    @test compare_pairlists(nlist_legacy, nlist_sort)
end

"""
    test_all_pbc_legacy_vs_sortbased(; N=50, density=0.05, cutoff_frac=0.33)

Run legacy vs sort-based comparison for all 8 PBC combinations.
"""
function test_all_pbc_legacy_vs_sortbased(; N=50, density=0.05, cutoff_frac=0.33)
    X, C, L = rand_config(N; density=density)
    cutoff = L * cutoff_frac

    for pbc in all_pbc_cases()
        nlist_legacy = PairList(X, cutoff, C, pbc; int_type=Int32)
        clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
        nlist_sort = materialize_pairlist(clist)

        @test npairs(nlist_legacy) == npairs(nlist_sort)
        @test compare_pairlists(nlist_legacy, nlist_sort)
    end
end

"""
    test_cpu_vs_gpu(to_gpu; N=100, density=0.05, cutoff_frac=0.33, pbc=SVec(true,true,true), T=Float64)

Test that GPU implementation matches CPU implementation.
`to_gpu` is a function that converts arrays to GPU (e.g., CuArray).
"""
function test_cpu_vs_gpu(to_gpu; N=100, density=0.05, cutoff_frac=0.33,
                          pbc=SVec(true, true, true), T=Float64)
    X, C, L = rand_config(N; density=density, T=T)
    cutoff = L * cutoff_frac

    # CPU version
    clist_cpu = build_cell_list(X, cutoff, C, pbc; backend=CPU())
    nlist_cpu = materialize_pairlist(clist_cpu)

    # GPU version
    X_gpu = to_gpu(X)
    clist_gpu = build_cell_list(X_gpu, cutoff, C, pbc)
    nlist_gpu = materialize_pairlist(clist_gpu)

    @test npairs(nlist_cpu) == npairs(nlist_gpu)
    @test compare_cpu_gpu_full(nlist_cpu, nlist_gpu)
end

"""
    test_all_pbc_cpu_vs_gpu(to_gpu; N=50, density=0.05, cutoff_frac=0.33, T=Float64)

Run CPU vs GPU comparison for all 8 PBC combinations.
"""
function test_all_pbc_cpu_vs_gpu(to_gpu; N=50, density=0.05, cutoff_frac=0.33, T=Float64)
    X, C, L = rand_config(N; density=density, T=T)
    cutoff = L * cutoff_frac

    for pbc in all_pbc_cases()
        clist_cpu = build_cell_list(X, cutoff, C, pbc; backend=CPU())
        nlist_cpu = materialize_pairlist(clist_cpu)

        X_gpu = to_gpu(X)
        clist_gpu = build_cell_list(X_gpu, cutoff, C, pbc)
        nlist_gpu = materialize_pairlist(clist_gpu)

        @test npairs(nlist_cpu) == npairs(nlist_gpu)
        @test compare_cpu_gpu_full(nlist_cpu, nlist_gpu)
    end
end

"""
    test_triclinic_cell(test_fn; N=80, cutoff=3.0)

Test with a non-cubic (triclinic) cell.
test_fn(X, C, cutoff, pbc) should run the actual test.
"""
function test_triclinic_cell(test_fn; N=80, cutoff=3.0, pbc=SVec(true, true, true),
                             T=Float64)
    C = rand_non_cubic_cell(T)
    X = [C' * SVec{T}(rand(T), rand(T), rand(T)) for _ in 1:N]
    test_fn(X, C, cutoff, pbc)
end

"""
    test_sizes(test_fn; sizes=[50, 100, 200, 500], density=0.05, cutoff_frac=0.25)

Run test_fn for multiple system sizes.
test_fn(X, C, L, cutoff) should run the actual test.
"""
function test_sizes(test_fn; sizes=[50, 100, 200, 500], density=0.05, cutoff_frac=0.25)
    for N in sizes
        X, C, L = rand_config(N; density=density)
        cutoff = L * cutoff_frac
        test_fn(X, C, L, cutoff)
    end
end

"""
    run_standard_correctness_suite(; gpu_array_fn=nothing)

Run the standard correctness test suite. If gpu_array_fn is provided,
also runs GPU tests.
"""
function run_standard_correctness_suite(; gpu_array_fn=nothing)
    @testset "Random configs" begin
        for _ in 1:5
            test_legacy_vs_sortbased(N=rand(50:200))
        end
    end

    @testset "All PBC combinations" begin
        test_all_pbc_legacy_vs_sortbased()
    end

    @testset "Triclinic cell" begin
        test_triclinic_cell() do X, C, cutoff, pbc
            nlist_legacy = PairList(X, cutoff, C, pbc; int_type=Int32)
            clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
            nlist_sort = materialize_pairlist(clist)
            @test npairs(nlist_legacy) == npairs(nlist_sort)
            @test compare_pairlists(nlist_legacy, nlist_sort)
        end
    end

    @testset "Large systems" begin
        test_sizes(sizes=[500, 1000]) do X, C, L, cutoff
            pbc = SVec(true, true, true)
            nlist_legacy = PairList(X, cutoff, C, pbc; int_type=Int32)
            clist = build_cell_list(X, cutoff, C, pbc; backend=CPU())
            nlist_sort = materialize_pairlist(clist)
            @test npairs(nlist_legacy) == npairs(nlist_sort)
        end
    end

    if gpu_array_fn !== nothing
        @testset "GPU correctness" begin
            test_cpu_vs_gpu(gpu_array_fn)
        end

        @testset "GPU all PBC" begin
            test_all_pbc_cpu_vs_gpu(gpu_array_fn)
        end
    end
end

# ==================== Edge Case Test Runners ====================

"""
    test_edge_cases(build_fn)

Run standard edge case tests using the provided build function.
`build_fn(X, cutoff, C, pbc)` should return a PairList.
"""
function test_edge_cases(build_fn)
    C = cubic_cell(10.0)

    # Single atom - no pairs
    nlist = build_fn([SVec(5.0, 5.0, 5.0)], 3.0, C, FULL_PBC)
    @test npairs(nlist) == 0

    # Two atoms within cutoff
    nlist = build_fn([SVec(5.0, 5.0, 5.0), SVec(5.0, 5.0, 6.0)], 3.0, C, FULL_PBC)
    @test npairs(nlist) == 2

    # Two atoms outside cutoff (no PBC)
    nlist = build_fn([SVec(1.0, 1.0, 1.0), SVec(8.0, 8.0, 8.0)], 3.0, C, NO_PBC)
    @test npairs(nlist) == 0
end

"""
    test_edge_cases_cpu_vs_gpu(to_gpu; T=Float64)

Run edge case tests comparing CPU and GPU implementations.
`to_gpu` is a function that converts arrays to GPU (e.g., CuArray).
"""
function test_edge_cases_cpu_vs_gpu(to_gpu; T=Float64)
    C = cubic_cell(T(10))

    # Single atom
    X = [SVec{T}(5, 5, 5)]
    nlist_gpu = materialize_pairlist(build_cell_list(to_gpu(X), 3.0, C, FULL_PBC))
    @test npairs(nlist_gpu) == 0

    # Two atoms - compare CPU vs GPU
    X = [SVec{T}(5, 5, 5), SVec{T}(5, 5, 6)]
    nlist_cpu = materialize_pairlist(build_cell_list(X, 3.0, C, FULL_PBC; backend=CPU()))
    nlist_gpu = materialize_pairlist(build_cell_list(to_gpu(X), 3.0, C, FULL_PBC))
    @test npairs(nlist_cpu) == npairs(nlist_gpu)
end

"""
    test_large_systems(build_fn; sizes=[500, 1000, 2000])

Run tests with various system sizes.
`build_fn(X, cutoff, C, pbc)` should return a PairList.
"""
function test_large_systems(build_fn; sizes=[500, 1000, 2000])
    for N in sizes
        X, C, L = rand_config(N)
        cutoff = L * 0.25
        nlist = build_fn(X, cutoff, C, FULL_PBC)
        @test npairs(nlist) > 0
    end
end

"""
    test_large_systems_cpu_vs_gpu(to_gpu; sizes=[500, 1000, 2000], T=Float64)

Run large system tests comparing CPU and GPU.
"""
function test_large_systems_cpu_vs_gpu(to_gpu; sizes=[500, 1000, 2000], T=Float64)
    for N in sizes
        X, C, L = rand_config(N; T=T)
        cutoff = L * 0.25
        nlist_cpu = materialize_pairlist(build_cell_list(X, cutoff, C, FULL_PBC; backend=CPU()))
        nlist_gpu = materialize_pairlist(build_cell_list(to_gpu(X), cutoff, C, FULL_PBC))
        @test npairs(nlist_cpu) == npairs(nlist_gpu)
    end
end

# ==================== Lazy Iteration Utilities ====================

"""
    count_lazy_pairs(clist)

Count total pairs by iterating through a SortedCellList lazily.
Useful for verifying lazy vs materialized pair counts match.
"""
function count_lazy_pairs(clist)
    count = 0
    for i in 1:nsites(clist)
        for_each_neighbour(clist, i) do j, R, S
            count += 1
        end
    end
    return count
end

"""
    count_manual_neighbours(nlist)

Count neighbours per atom by manually iterating through pair list.
Returns a vector where `result[i]` is the number of neighbours for atom i.
"""
function count_manual_neighbours(nlist)
    counts = zeros(Int, nsites(nlist))
    for idx in 1:npairs(nlist)
        counts[nlist.i[idx]] += 1
    end
    return counts
end

"""
    verify_neighbours_consistency(nlist, clist)

Verify that PairList and SortedCellList return the same neighbours for each atom.
Returns true if consistent, false otherwise.
"""
function verify_neighbours_consistency(nlist, clist)
    for i in 1:nsites(nlist)
        j_pairlist, _ = neigs(nlist, i)
        j_celllist, _, _ = neighbours(clist, i)
        if sort(j_pairlist) != sort(j_celllist)
            return false
        end
    end
    return true
end
