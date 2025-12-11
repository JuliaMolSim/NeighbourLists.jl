# GPU Kernels for NeighbourLists.jl
# Uses KernelAbstractions.jl for portable GPU programming

using KernelAbstractions
using KernelAbstractions: @kernel, @index, @Const, synchronize
using Atomix: @atomic

export compute_cell_ids_kernel!, count_neighbours_kernel!, fill_pairs_kernel!

# ==================== GPU Helper Functions ====================

"""
Convert 3D cartesian index to linear index (1-based, column-major).
"""
@inline function _sub2ind_gpu(dims::NTuple{3,TI}, idx::NTuple{3,TI}) where TI
    return idx[1] + (idx[2] - one(TI)) * dims[1] +
           (idx[3] - one(TI)) * dims[1] * dims[2]
end

@inline function _sub2ind_gpu(dims::SVec{TI}, idx::SVec{TI}) where TI
    return _sub2ind_gpu(dims.data, idx.data)
end

"""
Wrap cell index for periodic boundary and compute shift.
Returns (wrapped_index, shift).
"""
@inline function _wrap_and_shift_gpu(i::TI, n::TI, pbc::Bool) where TI
    if !pbc
        # For non-periodic: clamp to bounds, no shift
        return clamp(i, one(TI), n), zero(TI)
    end

    # For periodic: wrap around and track shift
    shift = zero(TI)
    wrapped = i

    while wrapped <= zero(TI)
        wrapped += n
        shift -= one(TI)
    end
    while wrapped > n
        wrapped -= n
        shift += one(TI)
    end

    return wrapped, shift
end

"""
Convert position to cell index (GPU-compatible version).
"""
@inline function _position_to_cell_index_gpu(inv_cell::SMat{T}, x::SVec{T},
                                              ncells::SVec{TI}) where {T, TI}
    # Transform to fractional coordinates and scale by cell count
    frac = inv_cell' * x
    return SVec{TI}(
        floor(TI, frac[1] * ncells[1] + one(TI)),
        floor(TI, frac[2] * ncells[2] + one(TI)),
        floor(TI, frac[3] * ncells[3] + one(TI))
    )
end

"""
Apply wrap or truncation to a cell index vector (GPU-compatible).
"""
@inline function _bin_wrap_or_trunc_vec_gpu(ci::SVec{TI}, pbc::SVec{Bool},
                                            ncells::SVec{TI}) where TI
    c1 = pbc[1] ? _bin_wrap_gpu(ci[1], ncells[1]) : clamp(ci[1], one(TI), ncells[1])
    c2 = pbc[2] ? _bin_wrap_gpu(ci[2], ncells[2]) : clamp(ci[2], one(TI), ncells[2])
    c3 = pbc[3] ? _bin_wrap_gpu(ci[3], ncells[3]) : clamp(ci[3], one(TI), ncells[3])
    return SVec{TI}(c1, c2, c3)
end

@inline function _bin_wrap_gpu(i::TI, n::TI) where TI
    while i <= zero(TI); i += n; end
    while i > n; i -= n; end
    return i
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

        # Transform to fractional coordinates and scale by cell count
        frac = inv_cell' * xi
        c = SVec{TI}(
            floor(TI, frac[1] * ncells[1] + one(TI)),
            floor(TI, frac[2] * ncells[2] + one(TI)),
            floor(TI, frac[3] * ncells[3] + one(TI))
        )

        # Apply boundary conditions
        c = _bin_wrap_or_trunc_vec_gpu(c, pbc, ncells)

        # Convert to linear index
        ns = ncells.data
        cell_ids[i] = _sub2ind_gpu(ns, c.data)
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

    # Use if-block instead of early return (KernelAbstractions requirement)
    if i <= nat
        xi = X_orig[i]
        TI = eltype(ncells)
        count = zero(TI)

        # Get cell index of atom i
        ci = _position_to_cell_index_gpu(inv_cell, xi, ncells)
        ci = _bin_wrap_or_trunc_vec_gpu(ci, pbc, ncells)

        # Iterate over adjacent cells
        for dz in -nxyz[3]:nxyz[3]
            for dy in -nxyz[2]:nxyz[2]
                for dx in -nxyz[1]:nxyz[1]
                    # Compute neighbor cell index with wrapping
                    cj_x, shift_x = _wrap_and_shift_gpu(ci[1] + TI(dx), ncells[1], pbc[1])
                    cj_y, shift_y = _wrap_and_shift_gpu(ci[2] + TI(dy), ncells[2], pbc[2])
                    cj_z, shift_z = _wrap_and_shift_gpu(ci[3] + TI(dz), ncells[3], pbc[3])

                    # Check if inside domain (for non-periodic)
                    in_bounds = true
                    if !pbc[1] && (ci[1] + TI(dx) < one(TI) || ci[1] + TI(dx) > ncells[1])
                        in_bounds = false
                    end
                    if !pbc[2] && (ci[2] + TI(dy) < one(TI) || ci[2] + TI(dy) > ncells[2])
                        in_bounds = false
                    end
                    if !pbc[3] && (ci[3] + TI(dz) < one(TI) || ci[3] + TI(dz) > ncells[3])
                        in_bounds = false
                    end

                    if in_bounds
                        cell_shift = SVec{TI}(shift_x, shift_y, shift_z)
                        cj_linear = _sub2ind_gpu((ncells[1], ncells[2], ncells[3]),
                                                 (cj_x, cj_y, cj_z))

                        # Iterate over atoms in this cell
                        idx_start = cell_offsets[cj_linear]
                        idx_end = cell_offsets[cj_linear + one(TI)] - one(TI)

                        for idx in idx_start:idx_end
                            j_orig = perm[idx]

                            # Skip self-interaction (unless periodic image)
                            is_self = (i == j_orig && shift_x == zero(TI) &&
                                       shift_y == zero(TI) && shift_z == zero(TI))

                            if !is_self
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

    # Use if-block instead of early return (KernelAbstractions requirement)
    if i <= nat
        xi = X_orig[i]
        TI = eltype(ncells)

        # Where to write pairs for atom i
        write_idx = offsets[i]

        # Get cell index of atom i
        ci = _position_to_cell_index_gpu(inv_cell, xi, ncells)
        ci = _bin_wrap_or_trunc_vec_gpu(ci, pbc, ncells)

        # Iterate over adjacent cells (same pattern as count kernel)
        for dz in -nxyz[3]:nxyz[3]
            for dy in -nxyz[2]:nxyz[2]
                for dx in -nxyz[1]:nxyz[1]
                    cj_x, shift_x = _wrap_and_shift_gpu(ci[1] + TI(dx), ncells[1], pbc[1])
                    cj_y, shift_y = _wrap_and_shift_gpu(ci[2] + TI(dy), ncells[2], pbc[2])
                    cj_z, shift_z = _wrap_and_shift_gpu(ci[3] + TI(dz), ncells[3], pbc[3])

                    # Check if inside domain (for non-periodic)
                    in_bounds = true
                    if !pbc[1] && (ci[1] + TI(dx) < one(TI) || ci[1] + TI(dx) > ncells[1])
                        in_bounds = false
                    end
                    if !pbc[2] && (ci[2] + TI(dy) < one(TI) || ci[2] + TI(dy) > ncells[2])
                        in_bounds = false
                    end
                    if !pbc[3] && (ci[3] + TI(dz) < one(TI) || ci[3] + TI(dz) > ncells[3])
                        in_bounds = false
                    end

                    if in_bounds
                        cell_shift = SVec{TI}(shift_x, shift_y, shift_z)
                        cj_linear = _sub2ind_gpu((ncells[1], ncells[2], ncells[3]),
                                                 (cj_x, cj_y, cj_z))

                        idx_start = cell_offsets[cj_linear]
                        idx_end = cell_offsets[cj_linear + one(TI)] - one(TI)

                        for idx in idx_start:idx_end
                            j_orig = perm[idx]

                            is_self = (i == j_orig && shift_x == zero(TI) &&
                                       shift_y == zero(TI) && shift_z == zero(TI))

                            if !is_self
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
"""
@kernel function map_pairs_symmetric_kernel!(out, f_half, @Const(i_arr), @Const(j_arr),
        @Const(X), @Const(S_arr), @Const(C))

    n = @index(Global, Linear)

    i, j = i_arr[n], j_arr[n]

    # Only process pairs with i < j (symmetric)
    if i < j
        # Compute displacement R
        R = X[j] - X[i] + C' * S_arr[n]

        # Compute value (passed as half since we add to both i and j)
        val = f_half

        # Atomic accumulation to both sites
        @atomic out[i] += val
        @atomic out[j] += val
    end
end


"""
Kernel for anti-symmetric pair mapping (force-like).
Adds +f to site j and -f to site i.
"""
@kernel function map_pairs_antisym_kernel!(out, f_val, @Const(i_arr), @Const(j_arr),
        @Const(X), @Const(S_arr), @Const(C))

    n = @index(Global, Linear)

    i, j = i_arr[n], j_arr[n]

    if i < j
        R = X[j] - X[i] + C' * S_arr[n]

        @atomic out[j] += f_val
        @atomic out[i] -= f_val
    end
end


# ==================== Extension Points ====================

# These functions are implemented in the CUDA extension

"""
Compute pair offsets from counts using GPU prefix sum.
Must be implemented in GPU extension (e.g., NeighbourListsCUDAExt).
"""
function compute_pair_offsets_gpu end

"""
GPU-accelerated materialize_pairlist.
Dispatched based on array backend.
"""
function _materialize_pairlist_gpu end
