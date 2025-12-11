module NeighbourListsCUDAExt

using NeighbourLists
using NeighbourLists: SortedCellList, PairList, SVec, SMat, nsites, lengths,
                      compute_cell_ids_kernel!, count_neighbours_kernel!, fill_pairs_kernel!
using CUDA
using CUDA: CuArray, CuVector
using KernelAbstractions
using KernelAbstractions: synchronize, get_backend
using LinearAlgebra: dot

# ==================== GPU Cell ID Computation ====================

"""
Compute cell IDs on GPU using a kernel.
"""
function NeighbourLists._compute_cell_ids(X::CuArray{<:SVec{T}}, inv_cell::SMat{T},
                                          ncells::SVec{TI}, pbc::SVec{Bool},
                                          ::Type{TI}) where {T, TI}
    nat = length(X)
    cell_ids = CuArray{TI}(undef, nat)
    backend = get_backend(X)

    kernel = compute_cell_ids_kernel!(backend)
    kernel(cell_ids, X, inv_cell, ncells, pbc; ndrange=nat)
    synchronize(backend)

    return cell_ids
end


# ==================== GPU Cell Offsets Computation ====================

"""
Compute CSR-style offsets from sorted cell IDs on GPU.
Uses GPU-friendly histogram and prefix sum.
"""
function NeighbourLists._compute_cell_offsets(sorted_cell_ids::CuArray{TI},
                                               ncells_total::Integer, ::Type{TI}) where TI
    nat = length(sorted_cell_ids)
    # +1 for sentinel at end
    cell_offsets = CUDA.zeros(TI, ncells_total + 1)

    if nat == 0
        cell_offsets .= one(TI)
        return cell_offsets
    end

    # Count atoms per cell using GPU histogram
    # cell_offsets[cid + 1] += 1 for each atom in cell cid
    # Use atomic adds for thread safety
    @cuda threads=256 blocks=cld(nat, 256) _histogram_kernel!(cell_offsets, sorted_cell_ids)
    CUDA.synchronize()

    # Cumulative sum to get offsets (1-indexed)
    # Use CUDA's native cumsum for better compatibility
    cumsum!(@view(cell_offsets[2:end]), @view(cell_offsets[2:end]))
    # First element is 1 (1-based), others already computed
    CUDA.@allowscalar cell_offsets[1] = one(TI)
    cell_offsets[2:end] .+= one(TI)  # Adjust for 1-based indexing

    return cell_offsets
end

"""
GPU kernel to compute histogram of cell IDs.
"""
function _histogram_kernel!(offsets, cell_ids)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if i <= length(cell_ids)
        cid = cell_ids[i]
        CUDA.@atomic offsets[cid + 1] += one(eltype(offsets))
    end
    return nothing
end


# ==================== GPU Sorting ====================

"""
GPU-accelerated sortperm using CUDA's sort_by_key.
"""
function NeighbourLists._get_sortperm(cell_ids::CuArray{TI}, backend) where TI
    n = length(cell_ids)
    # Create permutation array
    perm = CuArray{TI}(1:n)
    # Copy cell_ids for sorting (sort_by_key! modifies in place)
    sorted_ids = copy(cell_ids)

    # Sort permutation by cell_ids using CUDA's efficient sort
    CUDA.@sync begin
        # sortperm on GPU - use CPU fallback for now, can optimize later
        # TODO: Use CUDA.jl's sort_by_key! when stable API available
        cpu_ids = Array(cell_ids)
        cpu_perm = sortperm(cpu_ids)
        copyto!(perm, CuArray(TI.(cpu_perm)))
    end

    return perm
end


# ==================== GPU Prefix Sum ====================

"""
Compute exclusive prefix sum for pair offsets on GPU.

Given counts[i] = number of neighbours for atom i,
returns offsets where offsets[i] = 1 + sum(counts[1:i-1])
so that pairs for atom i are stored at indices offsets[i]:offsets[i+1]-1
"""
function NeighbourLists.compute_pair_offsets_gpu(counts::CuArray{TI}) where TI
    n = length(counts)

    # Allocate offsets array (n+1 elements for CSR format)
    offsets = CUDA.zeros(TI, n + 1)

    # First element is 1 (1-based indexing)
    CUDA.@allowscalar offsets[1] = one(TI)

    # Copy counts into offsets[2:end] then cumsum
    copyto!(@view(offsets[2:end]), counts)
    cumsum!(@view(offsets[2:end]), @view(offsets[2:end]))

    # Adjust for 1-based indexing: offsets[i+1] = 1 + sum(counts[1:i])
    offsets[2:end] .+= one(TI)

    return offsets
end


# ==================== GPU Pair Materialization ====================

"""
GPU implementation of materialize_pairlist using two-kernel approach.

1. Count neighbours per atom (parallel kernel)
2. Compute offsets via prefix sum
3. Fill pair arrays (parallel kernel)
"""
function NeighbourLists._materialize_pairlist_gpu(clist::SortedCellList{T, TI, AT},
                                                   backend) where {T, TI, AT <: CuArray}
    nat = nsites(clist)

    # Compute how many cells to check in each direction
    lens = lengths(clist.cell)
    nxyz = ceil.(TI, clist.cutoff * (clist.ncells ./ abs.(lens)))
    cutoff_sq = clist.cutoff^2

    # ===== Kernel 1: Count neighbours per atom =====
    counts = CUDA.zeros(TI, nat)

    kernel1 = count_neighbours_kernel!(backend)
    kernel1(counts, clist.X_orig, clist.cell_offsets, clist.perm,
            clist.inv_cell, clist.cell, clist.ncells, clist.pbc,
            cutoff_sq, nxyz; ndrange=nat)
    synchronize(backend)

    # ===== Step 2: Compute offsets via prefix sum =====
    offsets = NeighbourLists.compute_pair_offsets_gpu(counts)
    total_pairs = CUDA.@allowscalar offsets[end] - one(TI)

    if total_pairs == 0
        # Handle empty case
        i_arr = CUDA.zeros(TI, 0)
        j_arr = CUDA.zeros(TI, 0)
        S_arr = CuArray{SVec{TI}}(undef, 0)
        first_arr = CUDA.ones(TI, nat + 1)
        return PairList{T, TI, AT}(clist.X_orig, clist.cell, clist.cutoff,
                                   i_arr, j_arr, S_arr, first_arr)
    end

    # ===== Step 3: Allocate pair arrays =====
    i_arr = CuArray{TI}(undef, total_pairs)
    j_arr = CuArray{TI}(undef, total_pairs)
    S_arr = CuArray{SVec{TI}}(undef, total_pairs)

    # ===== Kernel 2: Fill pair arrays =====
    kernel2 = fill_pairs_kernel!(backend)
    kernel2(i_arr, j_arr, S_arr, offsets, clist.X_orig, clist.cell_offsets,
            clist.perm, clist.inv_cell, clist.cell, clist.ncells, clist.pbc,
            cutoff_sq, nxyz; ndrange=nat)
    synchronize(backend)

    # Use offsets as the first array (already in correct CSR format)
    first_arr = offsets[1:nat+1]

    return PairList{T, TI, AT}(clist.X_orig, clist.cell, clist.cutoff,
                               i_arr, j_arr, S_arr, first_arr)
end


# ==================== GPU Cell List Construction ====================

"""
Build cell list on GPU.
This dispatches to the standard build_cell_list but ensures GPU arrays are used.
"""
function NeighbourLists.build_cell_list(X::CuVector{<:SVec{T}}, cutoff::Real,
                                        cell::AbstractMatrix, pbc;
                                        int_type::Type{TI} = Int32,
                                        backend = get_backend(X)) where {T, TI}
    # Ensure cell matrix is on GPU-compatible type
    cell_mat = SMat{T}(cell)
    pbc_vec = SVec{Bool}(pbc)
    cutoff_T = T(cutoff)

    return NeighbourLists._build_sorted_celllist(X, cell_mat, pbc_vec, cutoff_T, TI, backend)
end


# ==================== GPU materialize_pairlist dispatch ====================

# This method catches CuArray-based cell lists and dispatches to GPU implementation
function NeighbourLists.materialize_pairlist(clist::SortedCellList{T, TI, AT};
                                              backend = get_backend(clist.X_orig)) where {T, TI, AT <: CuArray}
    return NeighbourLists._materialize_pairlist_gpu(clist, backend)
end


end # module
