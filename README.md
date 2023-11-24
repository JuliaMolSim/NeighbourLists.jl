# NeighbourLists.jl

A Julia port and restructuring of the neighbourlist implemented in
[matscipy](https://github.com/libAtoms/matscipy) (with the authors' permission).
Single-threaded, the Julia version is faster than matscipy for small systems,
probably due  to the overhead of dealing with Python, but on large systems it is
tends to be slower (up to around a factor 2 for 100,000 atoms). 

The package is can be used stand-alone, with
[JuLIP.jl](https://github.com/libAtoms/JuLIP.jl), or with [AtomsBase.jl](https://github.com/JuliaMolSim/AtomsBase.jl). 

## Getting Started

```Julia
Pkg.add("NeighbourLists")
using NeighbourLists
?PairList
```

Please also look at the tests on how to use this package. Or just open an issue and
ask.
