# NeighbourLists.jl

[![Build Status](https://travis-ci.org/libAtoms/NeighbourLists.jl.svg?branch=master)](https://travis-ci.org/libAtoms/NeighbourLists.jl)

[![Coverage Status](https://coveralls.io/repos/libAtoms/NeighbourLists.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/libAtoms/NeighbourLists.jl?branch=master)

[![codecov.io](http://codecov.io/github/libAtoms/NeighbourLists.jl/coverage.svg?branch=master)](http://codecov.io/github/libAtoms/NeighbourLists.jl?branch=master)


A Julia port and restructuring of the neighbourlist implemented in
[matscipy](https://github.com/libAtoms/matscipy) (with the authors' permission).
Single-threaded, the Julia version is faster than matscipy for small systems,
probably due  to the overhead of dealing with Python, but on large systems it is
tends to be slower (up to around a factor 2 for 100,000 atoms). However, the
Julia version is also multi-threaded with reasonably good scaling.

The package is intended to be used with
[JuLIP.jl](https://github.com/libAtoms/JuLIP.jl), but can be used as
stand-alone.

## Getting Started

```Julia
Pkg.add("NeighbourLists")
using NeighbourLists
?PairList
```

## TODO

* fix the performance problem, it should be possible to be closer to C
* 2D (structure in place, just needs a little time investment and testing)
* finish the iterators
