# GPU Kernels for NeighbourLists.jl
# Uses KernelAbstractions.jl for portable GPU programming

using KernelAbstractions
using KernelAbstractions: @kernel, @index, @Const, synchronize
using Atomix: @atomic

export compute_cell_ids_kernel!, count_neighbours_kernel!, fill_pairs_kernel!

# Note: Helper functions (_sub2ind_scalar, wrap_and_shift, position_to_cell_index_scalar,
# bin_wrap_or_trunc_vec, bin_wrap) are defined in cell_list.jl and work in GPU kernel contexts.

# ==================== Shared Kernel Helpers ====================

"""
Check if cell offset (dx, dy, dz) is within bounds for non-periodic boundaries.
"""
@inline function _is_cell_in_bounds(ci::SVec{TI}, dx::TI, dy::TI, dz::TI,
                                     ncells::SVec{TI}, pbc::SVec{Bool}) where TI
    (pbc[1] || (one(TI) <= ci[1] + dx <= ncells[1])) &&
    (pbc[2] || (one(TI) <= ci[2] + dy <= ncells[2])) &&
    (pbc[3] || (one(TI) <= ci[3] + dz <= ncells[3]))
end

"""
Check if this is a self-interaction (same atom without periodic shift).
"""
@inline function _is_self_interaction(i::Integer, j::Integer,
                                       shift_x::TI, shift_y::TI, shift_z::TI) where TI
    i == j && shift_x == zero(TI) && shift_y == zero(TI) && shift_z == zero(TI)
end

"""
Get the linear cell index and cell shift for a neighboring cell.
Returns (cell_linear_index, cell_shift_vector).
"""
@inline function _get_neighbor_cell(ci::SVec{TI}, dx::TI, dy::TI, dz::TI,
                                     ncells::SVec{TI}, pbc::SVec{Bool}) where TI
    cj_x, shift_x = wrap_and_shift(ci[1] + dx, ncells[1], pbc[1])
    cj_y, shift_y = wrap_and_shift(ci[2] + dy, ncells[2], pbc[2])
    cj_z, shift_z = wrap_and_shift(ci[3] + dz, ncells[3], pbc[3])
    cell_shift = SVec{TI}(shift_x, shift_y, shift_z)
    cj_linear = _sub2ind_scalar((ncells[1], ncells[2], ncells[3]), (cj_x, cj_y, cj_z))
    return cj_linear, cell_shift
end

# ==================== GPU Kernels ====================

"""
Kernel to compute cell ID for each atom position.
"""
@kernel function compute_cell_ids_kernel!(cell_ids, @Const(X), @Const(inv_cell),
        @Const(ncells), @Const(pbc))

    i = @index(Global, Linear)
    nat = length(X)

    if i <= nat
        xi = X[i]
        TI = eltype(ncells)

        # Compute cell index
        c = position_to_cell_index_scalar(inv_cell, xi, ncells)
        c = bin_wrap_or_trunc_vec(c, pbc, ncells)

        # Convert to linear index
        cell_ids[i] = _sub2ind_scalar(ncells, c)
    end
end


"""
Kernel to count neighbours for each atom.
Each thread handles one atom and counts its neighbours within cutoff.
"""
@kernel function count_neighbours_kernel!(counts, @Const(X_orig), @Const(cell_offsets),
        @Const(perm), @Const(inv_cell), @Const(cell_mat), @Const(ncells),
        @Const(pbc), cutoff_sq, @Const(nxyz))

    i = @index(Global, Linear)
    nat = length(X_orig)

    if i <= nat
        xi = X_orig[i]
        TI = eltype(ncells)
        count = zero(TI)

        # Get cell index of atom i
        ci = position_to_cell_index_scalar(inv_cell, xi, ncells)
        ci = bin_wrap_or_trunc_vec(ci, pbc, ncells)

        # Iterate over adjacent cells
        for dz in -nxyz[3]:nxyz[3]
            for dy in -nxyz[2]:nxyz[2]
                for dx in -nxyz[1]:nxyz[1]
                    if _is_cell_in_bounds(ci, TI(dx), TI(dy), TI(dz), ncells, pbc)
                        cj_linear, cell_shift = _get_neighbor_cell(ci, TI(dx), TI(dy), TI(dz), ncells, pbc)

                        # Iterate over atoms in this cell
                        idx_start = cell_offsets[cj_linear]
                        idx_end = cell_offsets[cj_linear + one(TI)] - one(TI)

                        for idx in idx_start:idx_end
                            j_orig = perm[idx]

                            if !_is_self_interaction(i, j_orig, cell_shift[1], cell_shift[2], cell_shift[3])
                                xj = X_orig[j_orig]
                                R = xj - xi + cell_mat' * cell_shift
                                r_sq = dot(R, R)

                                if r_sq < cutoff_sq
                                    count += one(TI)
                                end
                            end
                        end
                    end
                end
            end
        end

        counts[i] = count
    end
end


"""
Kernel to fill pair arrays using pre-computed offsets.
Each thread handles one atom and writes its pairs to the pre-allocated arrays.
"""
@kernel function fill_pairs_kernel!(i_arr, j_arr, S_arr, @Const(offsets),
        @Const(X_orig), @Const(cell_offsets), @Const(perm), @Const(inv_cell),
        @Const(cell_mat), @Const(ncells), @Const(pbc), cutoff_sq, @Const(nxyz))

    i = @index(Global, Linear)
    nat = length(X_orig)

    if i <= nat
        xi = X_orig[i]
        TI = eltype(ncells)

        # Where to write pairs for atom i
        write_idx = offsets[i]

        # Get cell index of atom i
        ci = position_to_cell_index_scalar(inv_cell, xi, ncells)
        ci = bin_wrap_or_trunc_vec(ci, pbc, ncells)

        # Iterate over adjacent cells (same pattern as count kernel)
        for dz in -nxyz[3]:nxyz[3]
            for dy in -nxyz[2]:nxyz[2]
                for dx in -nxyz[1]:nxyz[1]
                    if _is_cell_in_bounds(ci, TI(dx), TI(dy), TI(dz), ncells, pbc)
                        cj_linear, cell_shift = _get_neighbor_cell(ci, TI(dx), TI(dy), TI(dz), ncells, pbc)

                        idx_start = cell_offsets[cj_linear]
                        idx_end = cell_offsets[cj_linear + one(TI)] - one(TI)

                        for idx in idx_start:idx_end
                            j_orig = perm[idx]

                            if !_is_self_interaction(i, j_orig, cell_shift[1], cell_shift[2], cell_shift[3])
                                xj = X_orig[j_orig]
                                R = xj - xi + cell_mat' * cell_shift
                                r_sq = dot(R, R)

                                if r_sq < cutoff_sq
                                    i_arr[write_idx] = TI(i)
                                    j_arr[write_idx] = j_orig
                                    S_arr[write_idx] = cell_shift
                                    write_idx += one(TI)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


# ==================== GPU Map Operations ====================

"""
Kernel for symmetric pair mapping with atomic accumulation.
Each thread processes one pair (i,j) with i < j.
Calls f(i, j, R) and accumulates f/2 to both sites.
"""
@kernel function map_pairs_symmetric_kernel!(f, out, @Const(i_arr), @Const(j_arr),
        @Const(X), @Const(S_arr), @Const(C))

    n = @index(Global, Linear)

    i, j = i_arr[n], j_arr[n]

    # Only process pairs with i < j (symmetric)
    if i < j
        # Compute displacement R
        R = X[j] - X[i] + C' * S_arr[n]

        # Compute value and divide by 2 for symmetric accumulation
        val = f(i, j, R)
        val_half = val / 2

        # Atomic accumulation to both sites
        @atomic out[i] += val_half
        @atomic out[j] += val_half
    end
end


"""
Kernel for non-symmetric pair mapping with atomic accumulation.
Each thread processes one pair, accumulating to site i only.
"""
@kernel function map_pairs_nonsym_kernel!(f, out, @Const(i_arr), @Const(j_arr),
        @Const(X), @Const(S_arr), @Const(C))

    n = @index(Global, Linear)

    i, j = i_arr[n], j_arr[n]

    # Compute displacement R
    R = X[j] - X[i] + C' * S_arr[n]

    # Compute value and accumulate to site i
    val = f(i, j, R)

    @atomic out[i] += val
end


"""
Kernel for anti-symmetric pair mapping (force-like) - scalar version.
Adds +f to site j and -f to site i.
"""
@kernel function map_pairs_antisym_kernel!(f, out::AbstractVector{T}, @Const(i_arr), @Const(j_arr),
        @Const(X), @Const(S_arr), @Const(C)) where {T<:Real}

    n = @index(Global, Linear)

    i, j = i_arr[n], j_arr[n]

    if i < j
        R = X[j] - X[i] + C' * S_arr[n]

        val = f(i, j, R)

        @atomic out[j] += val
        @atomic out[i] -= val
    end
end

"""
Kernel for anti-symmetric pair mapping (force-like) - 3D vector version.
Uses a flat (3, N) array layout for GPU-friendly atomic operations.
Adds +f to site j and -f to site i.

Note: `out_flat` should be a (3, N) matrix where column i holds the 3D vector for atom i.
"""
@kernel function map_pairs_antisym_vec_kernel!(f, out_flat, @Const(i_arr), @Const(j_arr),
        @Const(X), @Const(S_arr), @Const(C))

    n = @index(Global, Linear)

    i, j = i_arr[n], j_arr[n]

    if i < j
        R = X[j] - X[i] + C' * S_arr[n]

        val = f(i, j, R)

        # Component-wise atomic operations on flat array
        @atomic out_flat[1, j] += val[1]
        @atomic out_flat[2, j] += val[2]
        @atomic out_flat[3, j] += val[3]
        @atomic out_flat[1, i] -= val[1]
        @atomic out_flat[2, i] -= val[2]
        @atomic out_flat[3, i] -= val[3]
    end
end


# Note: map_sites! is not implemented as a GPU kernel because the standard API
# f(Rs) expects a dynamically-sized vector of displacement vectors, which is
# incompatible with GPU execution. For GPU workflows, use map_pairs! or
# map_pairs_d! which operate on individual pairs and are fully GPU-accelerated.


# ==================== Portable GPU Operations ====================
# These use AcceleratedKernels.jl and Atomix.jl for multi-vendor GPU support
# (CUDA, ROCm, Metal, oneAPI)

"""
Kernel to compute histogram of cell IDs using atomic increments.
Works on any GPU backend via Atomix.jl.
"""
@kernel function _histogram_kernel!(offsets, @Const(cell_ids))
    i = @index(Global, Linear)
    if i <= length(cell_ids)
        cid = cell_ids[i]
        @atomic offsets[cid + 1] += one(eltype(offsets))
    end
end

"""
    compute_pair_offsets_gpu(counts::AbstractVector{TI}, backend) -> offsets

Compute CSR-style offsets from neighbour counts using portable GPU operations.
Uses AcceleratedKernels.accumulate! for prefix sum.

Given counts[i] = number of neighbours for atom i,
returns offsets where offsets[i] = 1 + sum(counts[1:i-1])
so that pairs for atom i are stored at indices offsets[i]:offsets[i+1]-1
"""
function compute_pair_offsets_gpu(counts::AbstractVector{TI}, backend) where TI
    n = length(counts)
    offsets = similar(counts, TI, n + 1)
    fill!(offsets, zero(TI))

    # Copy counts into offsets[2:end]
    copyto!(@view(offsets[2:end]), counts)

    # Prefix sum using AcceleratedKernels (portable across all GPU backends)
    AcceleratedKernels.accumulate!(+, @view(offsets[2:end]), @view(offsets[2:end]); init=zero(TI))

    # Adjust for 1-based indexing
    offsets .+= one(TI)

    return offsets
end

"""
    _compute_cell_ids_gpu(X, inv_cell, ncells, pbc, TI, backend) -> cell_ids

Compute cell IDs for all atoms on GPU using KernelAbstractions kernel.
"""
function _compute_cell_ids_gpu(X::AbstractVector{<:SVec{T}}, inv_cell::SMat{T},
                               ncells::SVec{TI}, pbc::SVec{Bool},
                               ::Type{TI}, backend) where {T, TI}
    nat = length(X)
    cell_ids = similar(X, TI, nat)

    kernel = compute_cell_ids_kernel!(backend)
    kernel(cell_ids, X, inv_cell, ncells, pbc; ndrange=nat)
    synchronize(backend)

    return cell_ids
end

"""
    _compute_cell_offsets_gpu(sorted_cell_ids, ncells_total, TI, backend) -> cell_offsets

Compute CSR-style cell offsets on GPU using portable histogram and prefix sum.
"""
function _compute_cell_offsets_gpu(sorted_cell_ids::AbstractVector{TI},
                                   ncells_total::Integer, ::Type{TI}, backend) where TI
    nat = length(sorted_cell_ids)
    cell_offsets = similar(sorted_cell_ids, TI, ncells_total + 1)
    fill!(cell_offsets, zero(TI))

    if nat == 0
        fill!(cell_offsets, one(TI))
        return cell_offsets
    end

    # Histogram using portable atomic increments via Atomix
    kernel = _histogram_kernel!(backend)
    kernel(cell_offsets, sorted_cell_ids; ndrange=nat)
    synchronize(backend)

    # Prefix sum using AcceleratedKernels
    AcceleratedKernels.accumulate!(+, @view(cell_offsets[2:end]), @view(cell_offsets[2:end]); init=zero(TI))

    # Adjust for 1-based indexing
    cell_offsets .+= one(TI)

    return cell_offsets
end

"""
    _materialize_pairlist_gpu(clist::SortedCellList, backend) -> PairList

GPU implementation of materialize_pairlist using two-kernel approach.
Portable across all GPU backends via KernelAbstractions + AcceleratedKernels.

1. Count neighbours per atom (parallel kernel)
2. Compute offsets via prefix sum (AcceleratedKernels)
3. Fill pair arrays (parallel kernel)
"""
function _materialize_pairlist_gpu(clist::SortedCellList{T, TI, AT, VI},
                                   backend) where {T, TI, AT, VI}
    nat = nsites(clist)

    if nat == 0
        i_arr = similar(clist.X_orig, TI, 0)
        j_arr = similar(clist.X_orig, TI, 0)
        S_arr = similar(clist.X_orig, SVec{TI}, 0)
        first_arr = similar(clist.X_orig, TI, 1)
        fill!(first_arr, one(TI))
        return PairList{T, TI, AT, typeof(i_arr), typeof(S_arr)}(
            clist.X_orig, clist.cell, clist.cutoff,
            i_arr, j_arr, S_arr, first_arr)
    end

    # Compute how many cells to check in each direction
    lens = lengths(clist.cell)
    nxyz = ceil.(TI, clist.cutoff * (clist.ncells ./ abs.(lens)))
    cutoff_sq = clist.cutoff^2

    # ===== Kernel 1: Count neighbours per atom =====
    counts = similar(clist.X_orig, TI, nat)
    fill!(counts, zero(TI))

    kernel1 = count_neighbours_kernel!(backend)
    kernel1(counts, clist.X_orig, clist.cell_offsets, clist.perm,
            clist.inv_cell, clist.cell, clist.ncells, clist.pbc,
            cutoff_sq, nxyz; ndrange=nat)
    synchronize(backend)

    # ===== Step 2: Compute offsets via prefix sum =====
    offsets = compute_pair_offsets_gpu(counts, backend)

    # Get total pairs (need scalar access)
    total_pairs = _scalar_getindex(offsets, length(offsets)) - one(TI)

    if total_pairs == 0
        i_arr = similar(clist.X_orig, TI, 0)
        j_arr = similar(clist.X_orig, TI, 0)
        S_arr = similar(clist.X_orig, SVec{TI}, 0)
        first_arr = similar(clist.X_orig, TI, nat + 1)
        fill!(first_arr, one(TI))
        return PairList{T, TI, AT, typeof(i_arr), typeof(S_arr)}(
            clist.X_orig, clist.cell, clist.cutoff,
            i_arr, j_arr, S_arr, first_arr)
    end

    # ===== Step 3: Allocate pair arrays =====
    i_arr = similar(clist.X_orig, TI, total_pairs)
    j_arr = similar(clist.X_orig, TI, total_pairs)
    S_arr = similar(clist.X_orig, SVec{TI}, total_pairs)

    # ===== Kernel 2: Fill pair arrays =====
    kernel2 = fill_pairs_kernel!(backend)
    kernel2(i_arr, j_arr, S_arr, offsets, clist.X_orig, clist.cell_offsets,
            clist.perm, clist.inv_cell, clist.cell, clist.ncells, clist.pbc,
            cutoff_sq, nxyz; ndrange=nat)
    synchronize(backend)

    # Use offsets as the first array (already in correct CSR format)
    first_arr = offsets[1:nat+1]

    return PairList{T, TI, AT, typeof(i_arr), typeof(S_arr)}(
        clist.X_orig, clist.cell, clist.cutoff,
        i_arr, j_arr, S_arr, first_arr)
end

# Helper to get scalar from GPU array (copies single element to CPU)
function _scalar_getindex(arr::AbstractVector, i::Integer)
    # For GPU arrays, copy the single element to CPU
    return Array(@view arr[i:i])[1]
end


# ==================== GPU Dispatch Methods ====================
# These are defined here (after cell_list.jl) so they can reference the GPU implementations

# GPU dispatch for _compute_cell_ids
function _compute_cell_ids(X::AbstractVector{SVec{T}}, inv_cell::SMat{T},
                           ncells::SVec{TI}, pbc::SVec{Bool}, ::Type{TI},
                           backend) where {T, TI}
    return _compute_cell_ids_gpu(X, inv_cell, ncells, pbc, TI, backend)
end

# GPU dispatch for _compute_cell_offsets
function _compute_cell_offsets(sorted_cell_ids::AbstractVector{TI},
                               ncells_total::Integer, ::Type{TI},
                               backend) where TI
    return _compute_cell_offsets_gpu(sorted_cell_ids, ncells_total, TI, backend)
end

# GPU dispatch for materialize_pairlist
function materialize_pairlist(clist::SortedCellList{T, TI, AT, VI};
                              backend = get_backend(clist.X_orig)) where {T, TI, AT, VI}
    return _materialize_pairlist_gpu(clist, backend)
end
