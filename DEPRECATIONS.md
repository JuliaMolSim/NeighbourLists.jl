# API Migration Guide

This document outlines the transition from the legacy linked-list algorithm to the unified sort-based implementation in NeighbourLists.jl.

## Version 0.6.x (Current)

The sort-based implementation is now the recommended API. The legacy linked-list implementation remains available as a reference implementation used internally for testing correctness.

### Legacy vs New API

| Legacy API | New API | Notes |
|------------|---------|-------|
| `PairList(X::Vector{SVec}, cutoff, cell, pbc)` | `neighbour_list(X, cutoff, cell, pbc)` | Sort-based, parallelizable |
| `CellList` struct | `SortedCellList` | Used internally by new API |
| `_celllist_` | `build_cell_list` | Internal function |
| `_pairlist_` | `materialize_pairlist` | Internal function |

> **Note:** The legacy implementation will be retained indefinitely as a reference implementation for validating correctness in tests. However, new code should use the unified `neighbour_list()` API.

### New Unified API

```julia
# High-level entry point (recommended)
nlist = neighbour_list(X, cutoff, cell, pbc)

# Lazy iteration (memory efficient)
clist = neighbour_list(X, cutoff, cell, pbc; lazy=true)
for_each_neighbour(clist, i) do j, R, S
    # process neighbour
end

# Unified accessors (work with both PairList and SortedCellList)
js, Rs, Ss = neighbours(nlist, i)
n = num_neighbours(nlist, i)
```

### AtomsBase Support

AtomsBase integration has been moved to a package extension. To use it:

```julia
using NeighbourLists
using AtomsBase, Unitful

# Extension loads automatically
nlist = PairList(system, 5.0u"Å")
```

## Why Use the New API?

### Benefits of Sort-Based Algorithm

1. **GPU support**: Works on CUDA, ROCm, Metal, and oneAPI via KernelAbstractions.jl
2. **Multi-threaded CPU**: Parallel construction and pair enumeration
3. **Memory efficiency**: Option for lazy iteration without materializing all pairs
4. **Consistent API**: Same code works on CPU and GPU

## Migration Guide

### Before (v0.5.x and earlier)

```julia
using NeighbourLists

# Legacy linked-list constructor
nlist = PairList(X, cutoff, cell, pbc)

# Access neighbours
j, R = neigs(nlist, i)
```

### After (v0.6.x+)

```julia
using NeighbourLists

# New unified API (recommended)
nlist = neighbour_list(X, cutoff, cell, pbc)

# Or explicitly with backend
nlist = neighbour_list(X, cutoff, cell, pbc; 
                       backend=NeighbourLists.CPU())

# GPU support
using CUDA
X_gpu = CuArray(X)
nlist_gpu = neighbour_list(X_gpu, cutoff, cell, pbc)

# Access neighbours (unchanged)
j, R = neigs(nlist, i)
# Or using unified accessor
j, R, S = neighbours(nlist, i)

# Lazy iteration (new, memory efficient)
clist = neighbour_list(X, cutoff, cell, pbc; lazy=true)
for_each_neighbour(clist, i) do j, R, S
    # process each neighbour
end
```

### AtomsBase Users

```julia
# Before: AtomsBase was always loaded
using NeighbourLists

# After: Load AtomsBase explicitly to enable extension
using NeighbourLists
using AtomsBase, Unitful

# Then use as before
nlist = PairList(system, 5.0u"Å")
clist = build_cell_list(system, 5.0u"Å")
```

## Questions?

If you have questions about migrating to the new API, please open an issue at:
https://github.com/JuliaMolSim/NeighbourLists.jl/issues
