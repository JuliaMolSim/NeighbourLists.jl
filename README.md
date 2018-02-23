# NeighbourList.jl

<!-- [![Build Status](https://travis-ci.org/cortner/NeighbourList.jl.svg?branch=master)](https://travis-ci.org/cortner/NeighbourList.jl)

[![Coverage Status](https://coveralls.io/repos/cortner/NeighbourList.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/cortner/NeighbourList.jl?branch=master)

[![codecov.io](http://codecov.io/github/cortner/NeighbourList.jl/coverage.svg?branch=master)](http://codecov.io/github/cortner/NeighbourList.jl?branch=master) -->

A Julia port of the neighbourlist implemented in
[matscipy](https://github.com/libAtoms/matscipy) (with the authors' permission).
The Julia version is faster than matscipy for small systems, probably due  to
the overhead of dealing with Python, but on large systems it is  about 30%
slower.

For now, see `test/runtests.jl` to look at how to use it.

### TODO

* multi-threading
* 2D
* iterators
* implement and compare against neighbourlist based on NearestNeighbours.jl
