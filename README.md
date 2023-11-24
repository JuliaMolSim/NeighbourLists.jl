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

### Usage via `AtomsBase.jl` 

```julia
using ASEconvert, NeighbourLists, Unitful
cu = ase.build.bulk("Cu") * pytuple((4, 2, 3))
sys = pyconvert(AbstractSystem, cu)
nlist = PairList(sys, 3.5u"Å")
neigs_1, Rs_1 = neigs(nlist, 1)
```

### Usage via `JuLIP.jl` 

```julia
using JuLIP 
at = bulk(:Cu) * (4, 2, 3)
nlist = neighbourlist(at, 3.5)
neigs_1, Rs_1 = neigs(nlist, 1)
``` 

Please also look at the tests on how to use this package. Or just open an issue and
ask.
