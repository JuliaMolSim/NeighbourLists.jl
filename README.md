# NeighbourLists.jl

A Julia port and restructuring of the neighbourlist implemented in
[matscipy](https://github.com/libAtoms/matscipy) (with the authors' permission).
Single-threaded, the Julia version is faster than matscipy for small systems,
probably due  to the overhead of dealing with Python, but on large systems it is
tends to be slower (up to around a factor 2 for 100,000 atoms). 

The package can be used stand-alone or with [AtomsBase.jl](https://github.com/JuliaMolSim/AtomsBase.jl).

## Getting Started

```Julia
Pkg.add("NeighbourLists")
using NeighbourLists
?PairList
```

### Usage via `AtomsBase.jl` 

```julia
using ASEconvert, NeighbourLists, Unitful
cu = ase.build.bulk("Cu") * pytuple((4, 2, 3))
sys = pyconvert(AbstractSystem, cu)
nlist = PairList(sys, 3.5u"Å")
neigs_1, Rs_1 = neigs(nlist, 1)
```

Please also look at the tests on how to use this package. Or just open an issue and
ask.

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
