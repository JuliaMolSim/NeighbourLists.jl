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

## GPU Support

GPU-accelerated neighbour list construction is available via CUDA.jl. When CUDA is loaded, the package automatically uses GPU kernels for systems stored on the GPU:

```julia
using LinearAlgebra, NeighbourLists, CUDA, StaticArrays

# Create positions on GPU
L = 10.0
X = CuArray([SVector{3,Float64}(L*rand(), L*rand(), L*rand()) for _ in 1:10000])
cell = SMatrix{3,3,Float64}(L*I)
pbc = SVector{3,Bool}(true, true, true)

# Build cell list and materialize pairs (runs on GPU)
clist = build_cell_list(X, 3.0, cell, pbc)
nlist = materialize_pairlist(clist)
```

The implementation uses a sort-based cell list with KernelAbstractions.jl for portable parallelism. On CPU this enables multi-threading; on GPU it uses native CUDA sorting and parallel kernels.

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
