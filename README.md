# NeighbourLists.jl

A Julia package for computing neighbour lists in molecular simulations. Originally a port of the neighbourlist from [matscipy](https://github.com/libAtoms/matscipy), now extended with multi-threaded CPU and portable GPU support.

The package can be used stand-alone or with [AtomsBase.jl](https://github.com/JuliaMolSim/AtomsBase.jl).

## Installation

```julia
using Pkg
Pkg.add("NeighbourLists")
```

## Two Implementations

The package provides two cell list implementations:

| Implementation | Algorithm | Parallelism | Use Case |
|---------------|-----------|-------------|----------|
| **Legacy** | Linked-list | Single-threaded | Small systems, compatibility |
| **Sort-based** | Sort by cell ID | Multi-threaded CPU, GPU | Large systems, performance |

Both produce identical results (validated in tests). The sort-based implementation is recommended for new code.

## Quick Start

### Using AtomsBase.jl

```julia
using AtomsBuilder, NeighbourLists, Unitful

sys = bulk(:Cu, cubic=true) * (4, 4, 4)
nlist = PairList(sys, 5.0u"Å")
j, R = neigs(nlist, 1)  # neighbours of atom 1
```

### Using raw arrays

```julia
using NeighbourLists, StaticArrays, LinearAlgebra

X = [SVector{3,Float64}(rand(3)...) for _ in 1:1000]  # positions
cell = SMatrix{3,3,Float64}(10.0*I)                    # 10x10x10 cell
pbc = SVector{3,Bool}(true, true, true)                # periodic

nlist = PairList(X, 3.0, cell, pbc)
j, R = neigs(nlist, 1)
```

## Sort-Based API (CPU and GPU)

The package provides a unified sort-based cell list that works on both CPU and GPU with the same API. The backend is automatically detected from the array type.

### CPU Example (Multi-threaded)

```julia
using NeighbourLists, StaticArrays, LinearAlgebra

# Create positions (CPU Vector)
L = 10.0
X = [SVector{3,Float64}(L*rand(), L*rand(), L*rand()) for _ in 1:10000]
cell = SMatrix{3,3,Float64}(L*I)
pbc = SVector{3,Bool}(true, true, true)

# Build cell list and materialize pairs
clist = build_cell_list(X, 3.0, cell, pbc)
nlist = materialize_pairlist(clist)

# Access neighbours of atom 1
j, R = neigs(nlist, 1)
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
clist = build_cell_list(X, 3.0, cell, pbc)
nlist = materialize_pairlist(clist)
```

**What's the same:** The `build_cell_list` and `materialize_pairlist` API is identical. Cell matrix, cutoff, and boundary conditions work the same way.

**What's different:** Only the array type changes (`Vector` vs `CuArray`/`ROCArray`/etc.). The backend is automatically detected - no need to specify it manually.

The implementation uses [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) for portable parallelism and [AcceleratedKernels.jl](https://github.com/JuliaGPU/AcceleratedKernels.jl) for portable sorting. On CPU this enables multi-threading; on GPU it runs native parallel kernels.

### Benchmarks

Benchmarks on NVIDIA RTX A4500 (cutoff = 5.0 Å, density = 0.05 atoms/Å³):

| Atoms | Pairs | CPU (1T) | CPU (16T) | GPU | GPU Speedup |
|------:|------:|---------:|----------:|--------:|------------:|
| 1,000 | 26k | 8 ms | 3.3 ms | 2.4 ms | 3.3x |
| 5,000 | 130k | 41 ms | 4.3 ms | 2.4 ms | 17x |
| 10,000 | 261k | 82 ms | 6.7 ms | 2.6 ms | 32x |
| 50,000 | 1.3M | 430 ms | 26 ms | 4.4 ms | 98x |
| 100,000 | 2.6M | 880 ms | 37 ms | 7.4 ms | 119x |

GPU throughput: ~360 million pairs/second for large systems.

*Note: GPU speedup is relative to CPU (1T). Run `scripts/benchmark.jl` with `-t 1` and `-t 16` to reproduce.*


### Acknowledgements

- Original inspiration from [matscipy](https://github.com/libAtoms/matscipy) neighbourlist written by Lars Pastewka
- Linked-list approach was implemented by Christoph Ortner
- Sort-based approach idea proposed by Teemu Järvinen and Timon Gutleb, and implemented by James Kermode
