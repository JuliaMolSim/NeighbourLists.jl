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
X = CuArray([SVector{3,Float64}(rand(3)...) for _ in 1:10000])
cell = SMatrix{3,3,Float64}(10I)
pbc = SVector{3,Bool}(true, true, true)

# Build cell list and materialize pairs (runs on GPU)
clist = build_cell_list(X, 3.0, cell, pbc)
nlist = materialize_pairlist(clist)
```

The implementation uses a sort-based cell list with KernelAbstractions.jl for portable parallelism. On CPU this enables multi-threading; on GPU it uses native CUDA sorting and parallel kernels.

### Benchmarks

Benchmarks on NVIDIA RTX A4500, 16 CPU threads (cutoff = 5.0 Å, density = 0.05 atoms/Å³):

| Atoms | Pairs | CPU (16T) | GPU | Speedup |
|------:|------:|----------:|--------:|--------:|
| 1,000 | 26k | 3.3 ms | 2.4 ms | 1.4x |
| 5,000 | 130k | 4.3 ms | 2.4 ms | 1.8x |
| 10,000 | 261k | 6.7 ms | 2.6 ms | 2.6x |
| 50,000 | 1.3M | 26 ms | 4.4 ms | 5.8x |
| 100,000 | 2.6M | 37 ms | 7.4 ms | 5.0x |

GPU throughput: ~360 million pairs/second for large systems.
