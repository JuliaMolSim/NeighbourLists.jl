# Deprecation Timeline

This document outlines the deprecation plan for NeighbourLists.jl as it transitions from the legacy linked-list algorithm to the unified sort-based implementation.

## Version 0.6.x (Current)

The legacy linked-list implementation remains functional but will emit deprecation warnings. Users are encouraged to migrate to the new API.

### Deprecated Functions and Types

| Deprecated | Replacement | Notes |
|------------|-------------|-------|
| `PairList(X::Vector{SVec}, cutoff, cell, pbc)` (linked-list) | `neighbour_list(X, cutoff, cell, pbc)` or `PairList(X::AbstractVector, cutoff, cell, pbc; backend=CPU())` | Sort-based, parallelizable |
| `CellList` struct | `SortedCellList` | Used internally by new API |
| `_celllist_` | `build_cell_list` | Internal function |
| `_pairlist_` | `materialize_pairlist` | Internal function |

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

## Version 0.7.0 (Planned)

The legacy linked-list implementation will be removed. This includes:

- Removal of `CellList` struct
- Removal of `_celllist_` and `_pairlist_` internal functions
- Removal of legacy `PairList(X::Vector{SVec}, ...)` constructor
- All code will use the sort-based algorithm exclusively

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
nlist = neighbour_list(X, cutoff, cell, pbc; backend=CPU())

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

## Suppressing Deprecation Warnings

If you need to suppress warnings during migration:

```julia
using Logging
with_logger(SimpleLogger(stderr, Logging.Error)) do
    # code using deprecated API
end
```

Or to suppress all deprecation warnings globally (not recommended for production):

```julia
Base.depwarn_mode[] = :none  # Suppress warnings
Base.depwarn_mode[] = :warn  # Re-enable warnings (default)
```

## Questions?

If you have questions about migrating to the new API, please open an issue at:
https://github.com/JuliaMolSim/NeighbourLists.jl/issues
