# NeighbourLists.jl

<!-- [![Build Status](https://travis-ci.org/cortner/NeighbourLists.jl.svg?branch=master)](https://travis-ci.org/cortner/NeighbourLists.jl)

[![Coverage Status](https://coveralls.io/repos/cortner/NeighbourLists.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/cortner/NeighbourLists.jl?branch=master)

[![codecov.io](http://codecov.io/github/cortner/NeighbourLists.jl/coverage.svg?branch=master)](http://codecov.io/github/cortner/NeighbourLists.jl?branch=master) -->

A Julia port of the neighbourlist implemented in
[matscipy](https://github.com/libAtoms/matscipy) (with the authors' permission).
The Julia version is faster than matscipy for small systems, probably due  to
the overhead of dealing with Python, but on large systems it is  about 30%
slower.

For now, see `test/runtests.jl` to look at how to use it.

### TODO

* 2D
* iterators
* multi-threading for neighbourlist construction and for iterators
* implement and compare against neighbourlist based on NearestNeighbours.jl
