# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NeighbourLists.jl is a Julia package for computing neighbour lists in molecular simulations. Originally a port of the matscipy neighbourlist, it now supports multi-threaded CPU and GPU (CUDA) execution via KernelAbstractions.jl. The package can be used standalone, with JuLIP.jl, or with AtomsBase.jl.

## Build and Test Commands

```bash
# Run all tests
julia --project -e 'using Pkg; Pkg.test()'

# Run tests in REPL
julia --project
]test

# Run a specific test file
julia --project test/test_celllist.jl

# Start Julia with the package loaded
julia --project -e 'using NeighbourLists'
```

## Architecture

### Core Data Structures (`src/types.jl`)
- `PairList{T,TI,AT}`: Main neighbour list storing pairs (i,j) with cell shifts S. Uses CSR-like `first` array for efficient per-atom neighbour lookup. `AT` parameterizes array type for GPU support.
- `SortedCellList{T,TI,AT}`: GPU-friendly cell list using sort-based construction. Atoms sorted by cell ID with CSR-style `cell_offsets` for O(1) cell access.
- `CellList{T,TI}`: Legacy linked-list cell structure (CPU only).
- `SVec{T}` and `SMat{T}`: Type aliases for 3D static vectors/matrices.

### New API (v0.6+): Sort-Based Cell List (`src/cell_list.jl`)
- `build_cell_list(X, cutoff, cell, pbc; backend=CPU())`: Build GPU-compatible cell list
- `materialize_pairlist(clist)`: Convert cell list to PairList (two-pass: count then fill)
- `for_each_neighbour(f, clist, i)`: Lazy iteration - call `f(j, R, S)` for each neighbour
- `get_neighbours(clist, i)`: Get all neighbours as vectors

### Legacy API (`src/cell_list.jl`)
- `_celllist_`: Bins atoms into cells using linked lists (sequential)
- `_pairlist_`: Constructs pairs by traversing linked lists
- `_fix_cell_`: Handles atoms outside cell boundaries

### Iteration Interfaces (`src/iterators.jl`)
- `pairs(nlist)`: Iterator over all (i, j, R) pairs
- `sites(nlist)`: Iterator over sites returning (i, j_neighbours, R_neighbours)

### MapReduce Operations (`src/mapreduce.jl`)
- `map_sites!(f, out, nlist)`: Apply function to each site's neighbourhood
- `map_pairs!(f, out, nlist)`: Apply function to pairs with thread-local accumulation
- `map_pairs_d!(f, out, nlist)`: Anti-symmetric variant for force-like quantities
- Legacy: `maptosites!`, `maptosites_d!` for iterator-based API

### AtomsBase Integration (`src/atoms_base.jl`)
- `PairList(system::AbstractSystem, cutoff; backend=CPU())`: Constructor for AtomsBase systems
- `build_cell_list(system::AbstractSystem, cutoff; backend=CPU())`: Lazy cell list from AtomsBase

### CUDA Extension (`ext/NeighbourListsCUDAExt.jl`)
- Loaded automatically when CUDA.jl is available
- Provides GPU-specific sorting

## Key Implementation Notes

- **Sort-based construction**: Atoms sorted by cell ID enables coalesced GPU memory access
- **Lazy vs materialized**: Use `build_cell_list` + `for_each_neighbour` for memory efficiency, or `materialize_pairlist` for batched operations
- **Backend selection**: Pass `backend=CPU()` or GPU backend; auto-detected from array type
- **3D-specific**: Code assumes 3D in `_sub2ind`, `CartesianIndices`, `lengths`
- **Integer type `TI`**: Parameterized for large systems (Int32 default, can use Int64/Int128)
- **Cell shifts `S`**: Store periodic image information as integer vectors
