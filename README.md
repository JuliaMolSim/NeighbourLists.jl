# NeighbourLists.jl

A Julia package for computing neighbour lists in molecular simulations. Originally a port of the neighbourlist from [matscipy](https://github.com/libAtoms/matscipy), now extended with multi-threaded CPU and portable GPU support.

The package can be used stand-alone or with [AtomsBase.jl](https://github.com/JuliaMolSim/AtomsBase.jl).

## Installation

```julia
using Pkg
Pkg.add("NeighbourLists")
```

## Unified API (Recommended)

The `neighbour_list()` function provides a unified entry point that works on both CPU and GPU with the same API. The backend is automatically detected from the array type.

> **Deprecation Notice:** The legacy `PairList` constructor using linked-list algorithm is deprecated and will be removed in v0.7.0. Use `neighbour_list()` instead. See [DEPRECATIONS.md](DEPRECATIONS.md) for migration details.

### CPU Example (Multi-threaded)

```julia
using NeighbourLists, StaticArrays, LinearAlgebra

# Create positions (CPU Vector)
L = 10.0
X = [SVector{3,Float64}(L*rand(), L*rand(), L*rand()) for _ in 1:10000]
cell = SMatrix{3,3,Float64}(L*I)
pbc = SVector{3,Bool}(true, true, true)

# Build neighbour list (uses sort-based algorithm with multi-threading)
nlist = neighbour_list(X, 3.0, cell, pbc)

# Access neighbours of atom 1
j, R = neighbours(nlist, 1)
```

### GPU Example (CUDA, ROCm, Metal, oneAPI)

```julia
using NeighbourLists, StaticArrays, LinearAlgebra
using CUDA  # or AMDGPU, Metal, oneAPI

# Create positions on GPU (only difference: use CuArray)
L = 10.0
X = CuArray([SVector{3,Float64}(L*rand(), L*rand(), L*rand()) for _ in 1:10000])
cell = SMatrix{3,3,Float64}(L*I)
pbc = SVector{3,Bool}(true, true, true)

# Same API - backend auto-detected from array type
nlist = neighbour_list(X, 3.0, cell, pbc)
```

**What's the same:** The `neighbour_list()` API is identical on CPU and GPU. Cell matrix, cutoff, and boundary conditions work the same way.

**What's different:** Only the array type changes (`Vector` vs `CuArray`/`ROCArray`/etc.). The backend is automatically detected - no need to specify it manually.

### Lazy Mode (Memory Efficient)

For large systems where materializing all pairs is memory-intensive, use lazy mode:

```julia
# Returns a SortedCellList instead of materializing all pairs
clist = neighbour_list(X, 3.0, cell, pbc; lazy=true)

# Iterate without storing all pairs in memory
for i in 1:nsites(clist)
    for_each_neighbour(clist, i) do j, R, S
        # process neighbour j with distance vector R and shift S
    end
end
```

### AtomsBase.jl Integration

```julia
using AtomsBuilder, NeighbourLists, Unitful

sys = bulk(:Cu, cubic=true) * (4, 4, 4)
nlist = PairList(sys, 5.0u"Å")
j, R = neighbours(nlist, 1)  # neighbours of atom 1
```

The implementation uses [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) for portable parallelism and [AcceleratedKernels.jl](https://github.com/JuliaGPU/AcceleratedKernels.jl) for portable sorting. On CPU this enables multi-threading; on GPU it runs native parallel kernels.

## Two Implementations

The package provides two cell list implementations:

| Implementation | Algorithm | Parallelism | Status |
|---------------|-----------|-------------|--------|
| **Sort-based** | Sort by cell ID | Multi-threaded CPU, GPU | Recommended |
| **Legacy** | Linked-list | Single-threaded | Deprecated in v0.6, removed in v0.7 |

Both produce identical results (validated in tests). The sort-based implementation accessed via `neighbour_list()` is recommended for all new code.

## Deprecation Notice

The legacy linked-list implementation (`CellList`, `_celllist_`, `_pairlist_`) is deprecated and will be removed in v0.7.0.

**Migration summary:**
- Replace `PairList(X, cutoff, cell, pbc)` with `neighbour_list(X, cutoff, cell, pbc)`
- Replace `neigs(nlist, i)` with `neighbours(nlist, i)` (both still work)
- For memory-efficient iteration, use `neighbour_list(...; lazy=true)` with `for_each_neighbour`

See [DEPRECATIONS.md](DEPRECATIONS.md) for the complete migration guide.

## Benchmarks

Benchmarks on NVIDIA RTX A4500 (cutoff = 5.0 Å, density = 0.05 atoms/Å³):

| Atoms | Pairs | Legacy | CPU (1T) | CPU (8T) | GPU | Speedup |
|------:|------:|-------:|---------:|---------:|--------:|--------:|
| 1,000 | 26k | 8 ms | 3.6 ms | 3.4 ms | 2.3 ms | 3.5x |
| 5,000 | 131k | 38 ms | 17 ms | 3.9 ms | 2.2 ms | 17x |
| 10,000 | 262k | 84 ms | 35 ms | 7.8 ms | 2.4 ms | 36x |
| 50,000 | 1.3M | 516 ms | 201 ms | 31 ms | 4.2 ms | 124x |
| 100,000 | 2.6M | 1.1 s | 400 ms | 62 ms | 6.9 ms | 160x |

GPU throughput: ~370 million pairs/second for large systems.

*Note: Speedup is GPU vs Legacy. Run `julia --project -t N scripts/benchmark.jl` to reproduce.*


## Acknowledgements

- Original inspiration from [matscipy](https://github.com/libAtoms/matscipy) neighbourlist written by Lars Pastewka
- Linked-list approach was implemented by Christoph Ortner
- Sort-based approach idea proposed by Teemu Järvinen and Timon Gutleb, and implemented by James Kermode
