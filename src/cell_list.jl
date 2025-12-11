using Base.Threads, LinearAlgebra

export npairs, nsites, maxneigs, max_neighbours, neigs, neighbours, neigs!
export build_cell_list, materialize_pairlist, for_each_neighbour, get_neighbours, count_neighbours

PairList(X::Vector{SVec{T}}, cutoff::AbstractFloat, cell::AbstractMatrix, pbc;
            int_type::Type = Int32, fixcell = true) where {T} =
   _pairlist_(X, SMat{T}(cell), SVec{Bool}(pbc), T(cutoff), int_type,
              fixcell)

PairList(X::Matrix{T}, args...; kwargs...) where {T} =
   PairList(reinterpret(SVec{T}, X, (size(X,2),)), args...; kwargs...)

npairs(nlist::PairList) = length(nlist.i)
nsites(nlist::PairList) = length(nlist.first) - 1
cutoff(nlist::PairList) = nlist.cutoff


# ====================== Cell Index Algebra =====================

"Map i back to the interval [0,n) by shifting by integer multiples of n"
@inline function bin_wrap(i::Integer, n::Integer)
   while i <= 0;  i += n; end
   while i > n;  i -= n; end
   return i;
end

"apply bin_wrap only if periodic bdry"
@inline bin_wrap(i::Integer, pbc::Bool, n::Integer) = pbc ? bin_wrap(i, n) : i

"Map i back to the interval [0,n) by assigning edge value if outside interval"
@inline function bin_trunc(i::TI, n::TI) where {TI <: Integer}
   if i <= 0;     i = TI(1)
   elseif i > n;  i = n
   end
   return i
end

"apply bin_trunc only if open bdry"
@inline bin_trunc(i::TI, pbc::Bool, n::TI) where {TI <: Integer} =
      pbc ? TI(i) : bin_trunc(i, n)

"applies bin_trunc to open bdry and bin_wrap to periodic bdry"
@inline bin_wrap_or_trunc(i::Integer, pbc::Integer, n::Integer) =
      pbc ? bin_wrap(i, n) : bin_trunc(i, n)

"Map particle position to a (cartesian) cell index"
@inline position_to_cell_index(inv_cell::SMat{T}, x::SVec{T}, ns::SVec{TI}
         ) where {T, TI <: Integer} =
   floor.(TI, ((inv_cell' * x) .* ns .+ 1))


# ------------ The next two functions are the only dimension-dependent
#              parts of the code!

# an extension of sub2ind for the case when i is a vector (cartesian index)
# @inline Base.sub2ind(dims::NTuple{3,TI}, i::SVec{TI}) where {TI <: Integer} =
#    sub2ind(dims, i[1], i[2], i[3])
# WARNING: this smells like a performance regression!
# WARNING: `LinearIndices` always returnsstores Int, not TI!!!
@inline _sub2ind(dims::NTuple{3,TI}, i::SVec{TI}) where {TI <: Integer} =
   TI( (LinearIndices(dims))[i[1], i[2], i[3]] )

lengths(C::SMat{T}) where {T} =
   det(C) ./ SVec{T}(norm(C[2,:]×C[3,:]), norm(C[3,:]×C[1,:]), norm(C[1,:]×C[2,:]))


# --------------------------------------------------------------------------

function analyze_cell(cell, cutoff, TI)
   @assert TI <: Integer
   # check the cell volume (allow only 3D volumes!)
   volume = abs(det(cell))
   if volume < 1e-12
      @warn("zero volume detected - proceed at your own risk")
      @show cell
      @show volume
      @show cutoff
   end
   # precompute inverse of cell matrix for coordiate transformation
   inv_cell = inv(cell)
   # Compute distance of cell faces
   lens = abs.(lengths(cell))
   # Number of cells for cell subdivision
   _t = floor.(TI, lens / cutoff)
   ns_vec = max.(_t, one(TI))
   return inv_cell, ns_vec, lens
end

# multi-threading setup

function setup_mt(niter::TI, maxnt = MAX_THREADS[1]) where TI <: Integer
   nt = minimum([6, nthreads(), ceil(TI, niter / 20), maxnt])
   # nn = ceil.(TI, linspace(1, niter+1, nt+1))
   # range(start, stop=stop, length=length)
   nn = ceil.(TI, range(1, stop=niter+1, length=nt+1))
   return nt, nn
end

# ==================== CellList Core  ================


function _celllist_(X::Vector{<:SVec}, cell::SMat, pbc::SVec{Bool},
            cutoff::T, TI) where {T <: Real}
   @assert TI <: Integer
   # ----- analyze cell -----
   nat = length(X)
   inv_cell, ns_vec, lens = analyze_cell(cell, cutoff, TI)
   ns = ns_vec.data

   if prod(BigInt.(ns_vec)) > typemax(TI)
      error("""Ratio of simulation cell size to cutoff is very large.
               Are you using a cell with lots of vacuum? To fix this
               use a larger integer type (e.g. Int128), a
               larger cut-off, or a smaller simulation cell.""")
   end

   # data structure to store a linked list for each bin
   ncells = prod(ns_vec)
   seed = fill(TI(-1), ncells)
   last = Vector{TI}(undef, ncells)
   next = Vector{TI}(undef, nat)
   nats = zeros(TI, ncells)

   for i = 1:nat
      # Get cell index
      c = position_to_cell_index(inv_cell, X[i], ns_vec)
      # Periodic/non-periodic boundary conditions
      c = bin_wrap_or_trunc.(c, pbc, ns_vec)
      # linear cell index  # (+1 due to 1-based indexing)
      ci = _sub2ind(ns, c)

      # Put atom into appropriate bin (list of linked lists)
      if seed[ci] < 0   #  ci contains no atom yet
         next[i] = -1
         seed[ci] = i
         last[ci] = i
         nats[ci] += 1
      else
         next[i] = -1
         next[last[ci]] = i
         last[ci] = i
         nats[ci] += 1
      end
   end

   return CellList(X, cell, inv_cell, pbc, cutoff, seed, last, next, nats)
end


function _pairlist_(clist::CellList{T, TI}) where {T, TI}

   X, cell, pbc, cutoff, seed, last, next =
         clist.X, clist.cell, clist.pbc, clist.cutoff,
         clist.seed, clist.last, clist.next
   nat = length(X)
   inv_cell, ns_vec, lens = analyze_cell(cell, cutoff, TI)

   # guess how many neighbours per atom
   #    atoms in cell x (8 cells) * (ball / cube)
   max_nat_cell = maximum(clist.nats)
   nneigs_guess = ceil(TI, 1.5 * max_nat_cell * (π^2/4))

   # allocate arrays for the many threads
   # set number of threads
   nt, nn = setup_mt(nat)
   # allocate arrays
   first_t = Vector{TI}[ Vector{TI}()  for n = 1:nt ]    # i
   secnd_t = Vector{TI}[ Vector{TI}()  for n = 1:nt ]    # j
   shift_t = Vector{SVec{TI}}[ Vector{SVec{TI}}()  for n = 1:nt ]  # ~ X[i] - X[j]
   # give size hints
   sz = (nat ÷ nt + nt) * nneigs_guess
   for n = 1:nt
      sizehint!(first_t[n], sz)
      sizehint!(secnd_t[n], sz)
      sizehint!(shift_t[n], sz)
   end

   # We need the shape of the bin ( bins[:, i] = cell[i,:] / ns[i] )
   # bins = cell' ./ ns_vec
   bins = hcat( cell[1,:]/ns_vec[1], cell[2,:]/ns_vec[2], cell[3,:] / ns_vec[3] )

   # Find out over how many neighbor cells we need to loop (if the box is small)
   nxyz = ceil.(TI, cutoff * (ns_vec ./ lens))
   # cxyz = CartesianIndex(nxyz.data)
   # WARNING : 3D-specific hack; also potential performance regression
   # WARNING: `CartesianIndices` always stores Int, not TI!!!
   xyz_range = CartesianIndices((-nxyz[1]:nxyz[1],
                                 -nxyz[2]:nxyz[2],
                                 -nxyz[3]:nxyz[3]))

   # Loop over threads
   # @threads for it = 1:nt
   for it = 1:nt
      for i = nn[it]:(nn[it+1]-1)
         # current atom position
         _find_neighbours_!(i, clist, ns_vec, bins, xyz_range,
                      first_t[it], secnd_t[it], shift_t[it])
      end # for i = 1:nat
   end # @threads

   # piece them back together
   sz = sum( length(first_t[i]) for i = 1:nt )
   first = first_t[1];     sizehint!(first, sz)
   secnd = secnd_t[1];     sizehint!(secnd, sz)
   shift = shift_t[1]; sizehint!(shift, sz)
   for it = 2:nt
      append!(first, first_t[it])
      append!(secnd, secnd_t[it])
      append!(shift, shift_t[it])
   end

   # Build return tuple
   return first, secnd, shift
end



function _find_neighbours_!(i, clist, ns_vec::SVec{TI}, bins, xyz_range,
                            first, secnd, shift) where TI
   inv_cell, X, pbc = clist.inv_cell, clist.X, clist.pbc
   seed, last, next = clist.seed, clist.last, clist.next
   xi = X[i]

   # funnily testing with cutoff^2 actually makes a measurable difference
   # which suggests we are pretty close to the performance limit
   cutoff_sq = clist.cutoff^2

   # cell index (cartesian) of xi
   ci0 = position_to_cell_index(inv_cell, xi, ns_vec)

   # Truncate if non-periodic and outside of simulation domain
   # (here, we don't yet want to wrap the pbc as well)
   ci = bin_trunc.(ci0, pbc, ns_vec)
   # dxi is the position relative to the lower left corner of the bin
   dxi = xi - bins * (ci .- 1)

   # Apply periodic boundary conditions as well now
   ci = bin_wrap_or_trunc.(ci0, pbc, ns_vec)

   for ixyz in xyz_range  # Integer
      # convert cartesian index to SVector
      xyz = SVec{TI}(ixyz.I)
      # get the bin index
      cj = bin_wrap.(ci + xyz, pbc, ns_vec)
      # skip this bin if not inside the domain
      all(TI(1) .<= cj .<= ns_vec) || continue
      # linear cell index
      ncj = _sub2ind(ns_vec.data, cj)
      # Offset of the neighboring bins
      off = bins * xyz

      # Loop over all atoms in neighbouring bin (all potential
      # neighbours in the bin with linear index cj1)
      j = seed[ncj] # the first atom in the ncj cell
      while j > 0
         if i != j || any(xyz .!= 0)
            xj = X[j] # position of current neighbour

            # we need to find the cell index again, because this is
            # not really the cell index, but it could be outside
            # the domain -> i.e. this only makes a difference for pbc
            cj = position_to_cell_index(inv_cell, xj, ns_vec)
            cj = bin_trunc.(cj, pbc, ns_vec)

            # drj is position relative to lower left corner of the bin
            dxj = xj - bins * (cj .- 1)
            # Compute distance between atoms
            dx = dxj - dxi + off
            norm_dx_sq = dot(dx, dx)

            # append to the list
            if norm_dx_sq < cutoff_sq
               push!(first, i)
               push!(secnd, j)
               push!(shift, (ci0 - cj + xyz) .÷ ns_vec)
            end
         end  # if i != j || any(xyz .!= 0)

         # go to the next atom in the current cell
         j = next[j];
      end # while j > 0 (loop over atoms in current cell)
   end # loop over neighbouring bins
end


function _pairlist_(X::Vector{SVec{T}}, cell::SMat{T}, pbc::SVec{Bool},
            cutoff::T, TI, fixcell::Bool) where {T}
   @assert TI <: Integer
   # temporary (?) fix to make sure all atoms are within the cell
   if fixcell
      X, cell = _fix_cell_(X, cell, pbc)
   end

   clist = _celllist_(X, cell, pbc, cutoff, TI)
   i, j, S = _pairlist_(clist)

   first = get_first(i, length(X))
   sort_neigs!(j, (S,), first)

   return PairList{T, TI, Vector{SVec{T}}, Vector{TI}, Vector{SVec{TI}}}(X, cell, cutoff, i, j, S, first)
end



# ==================== Post-Processing =========================

"""
`get_first(i::Vector{TI}, nat::Integer) -> first::Vector{TI}
where TI <: Integer`

Assumes that `i` is sorted in ascending order.
For `n = 1, . . ., nat`, `first[n]` will be the index to the first
element of `i` with value `n`. Further, `first[nat+1]` will be
`length(i) + 1`.

If `first[n] == first[n+1]` then this means that `i` contains no element `n`.
"""
function get_first(i::Vector{TI}, nat::Integer = i[end]) where {TI}
   # compute the first index for each site
   first = Vector{TI}(undef, nat + 1)
   idx = 1
   n = 1
   while n <= nat && idx <= length(i)
      first[n] = idx
      while i[idx] == n
         idx += 1
         if idx > length(i)
            break
         end
      end
      n += 1
   end
   first[n:end] .= length(i)+1
   return first
end


"""
`sort_neigs!(j, arrays, first)`

sorts each sub-range of `j` corresponding to one site  in ascending order
and applies the same permutation to `r, R, S`.
"""
function sort_neigs!(j, arrays::Tuple, first)
   nat = length(first) - 1
   nt, nn = setup_mt(nat)
   @threads for it = 1:nt
      for n = nn[it]:(nn[it+1]-1)
         if first[n+1] > first[n] + 1
            rg = first[n]:first[n+1]-1
            jrg = j[rg]
            if !issorted(jrg)
               I = sortperm(j[rg])
               rg_perm = rg[I]
               j[rg] = j[rg_perm]
               for a in arrays
                  a[rg] .= a[rg_perm]
               end
            end
         end
      end
   end
end


"""
`_fix_cell_(X::Vector{SVec{T}}, C::SMat{T}, pbc)`

produces new `X`, `C` such that PBC are respected, but all positions are
inside the cell. This Potentially involves
 * wrapping atom positions in the pbc directions
 * shifting atom positions in the non-pbc directions
 * enlarging the cell
If either `X` or `C` are modified, then they will be copied.
"""
function _fix_cell_(X::Vector{SVec{T}}, C::SMat{T}, pbc) where {T}
   invCt = inv(C)'
   copy_X = false
   min_lam = @MVector zeros(3)
   max_lam = @MVector ones(3)
   for (ix, x) in enumerate(X)
      λ = Vector(invCt * x)
      for i = 1:3
         update_x = false
         if !(0.0 <= λ[i] < 1.0)
            if pbc[i]
               λ[i] = mod(λ[i], 1.0)
               update_x = true
            else
               min_lam[i] = min(min_lam[i], λ[i])
               max_lam[i] = max(max_lam[i], λ[i])
            end
         end
         if update_x
            if !copy_X; X = copy(X); copy_X = true end
            X[ix] = C' * SVec{T}(λ)
         end
      end
   end
   # check whether we need to adjust the non-PBC directions
   if (minimum(min_lam) < 0.0) || (maximum(max_lam) > 1.0)
      # shift vector:
      if !copy_X; X = copy(X); copy_X = true end
      t = - C' * min_lam
      for n = 1:length(X)
         X[n] += t
      end
      # the new min_lam is now zero and the new max_lam is max_lam - min_lam
      max_lam -= min_lam
      min_lam .= 0.0
      # because of round-off we add a little extra in the non-pbc directions 
      # it doesn't affect any energies etc but creates weirdness otherwise.
      for i = 1:3 
         if max_lam[i] > 1
            max_lam[i] *= 1.01
         end
      end
      # we need to multiply the cell by the correct lambda
      C = hcat( max_lam[1] * C[1,:], max_lam[2] * C[2,:], max_lam[3] * C[3,:])'
   end
   return X, C
end

"""
`maxneigs(nlist::PairList) -> Integer`

returns the maximum number of neighbours that any atom in the
neighbourlist has.
"""
maxneigs(nlist::PairList) = maximum( nneigs(nlist, n) for n = 1:nsites(nlist) )

# retire max_neigs
const max_neigs = maxneigs

"""
`nneigs(nlist::PairList, i0::Integer) -> Integer` :
return number of neighbours of particle with index `i0`.
(within cutoff radius + buffer)
"""
nneigs(nlist::PairList, i0::Integer) = nlist.first[i0+1]-nlist.first[i0]

function _getR(nlist, idx)
   i = nlist.i[idx]
   j = nlist.j[idx]
   return _getR(nlist.X[j] - nlist.X[i], nlist.S[idx], nlist.C)
end

_getR(dX::SVec, S::SVec, C::SMat) = dX + C' * S

"""
`neigs!(Rtemp, nlist, i) -> j, R`

For `nlist::PairList` this returns the interaction neighbourhood of
the atom indexed by `i`. E.g., in the standard loop approach one
would have
```
Rtemp = zeros(JVecF, maxneigs(nlist))
for (i, j, R) in sites(nlist)
   (j, R) == neigs(nlist, i) == neigs!(Rtemp, nlist, i)
end
```

`R` is a view into `Rtemp` with the correct length, while `j` is a view
into `nlist.j`.
"""
function neigs!(Rs::AbstractVector{<: SVec}, nlist::PairList, i0::Integer)
   n1, n2 = nlist.first[i0], nlist.first[i0+1]-1
   _grow_array!(Rs, n2-n1+1)
   J = (@view nlist.j[n1:n2])
   for n = 1:length(J)
      Rs[n] = _getR(nlist, n1+n-1)
   end
   return J, (@view Rs[1:length(J)])
end

function neigs!(Js::AbstractVector{<: SVec},
                Rs::AbstractVector{<: SVec},
                nlist::PairList, i0::Integer)
   _grow_array!(Rs, n2-n1+1)
   _grow_array!(Js, n2-n1+1)
   j, Rs = neigs!(Rs, nlist, i0)
   N = length(j)
   copyto!(Js, j)
   return (@view Js[1:length(j)]), Rs
end

function _grow_array!(A::Vector{T}, N) where {T}
   if length(A) < N
      append!(A, zeros(T, N-length(A)))
   end
   return A
end

"""
`neigss!(Rs, nlist, i0) -> j, R, S` : return neighbourhood as in
`neigs!` as well as the corresponding cell shifts.

(`R` is a view into `Rs` with the correct length)
"""
function neigss!(Rs::AbstractVector{<: SVec}, nlist::PairList, i0::Integer)
   n1, n2 = nlist.first[i0], nlist.first[i0+1]-1
   _grow_array!(Rs, n2-n1+1)
   J = (@view nlist.j[n1:n2])
   for n = 1:length(J)
      Rs[n] = _getR(nlist, n1+n-1)
   end
   return J, (@view Rs[1:length(J)]), (@view nlist.S[n1:n2])
end


neigs(nlist::PairList{T}, i0::Integer) where {T} =
      neigs!( zeros(SVec{T}, nneigs(nlist, i0)), nlist, i0 )

neigss(nlist::PairList{T}, i0::Integer) where {T} =
      neigss!( zeros(SVec{T}, nneigs(nlist, i0)), nlist, i0 )

"""
alias for `neigs`
"""
neighbours = neigs

"""
alias for `max_neigs`
"""
max_neighbours = maxneigs


# ==================== Sort-Based Cell List (GPU-compatible) ====================

"""
    build_cell_list(X, cutoff, cell, pbc; int_type=Int32, backend=CPU())

Build a GPU-compatible cell list using sort-based construction.

Returns a `SortedCellList` where atoms are sorted by their cell ID,
enabling efficient parallel traversal.

# Arguments
- `X`: Vector of positions (SVec{T} or similar)
- `cutoff`: Cutoff distance
- `cell`: 3×3 cell matrix (rows are lattice vectors)
- `pbc`: Periodic boundary conditions (3-tuple or SVec{Bool})
- `int_type`: Integer type for indices (default Int32)
- `backend`: KernelAbstractions backend (default CPU())
"""
function build_cell_list(X::AbstractVector{<:SVec{T}}, cutoff::Real,
                         cell::AbstractMatrix, pbc;
                         int_type::Type{TI} = Int32,
                         backend = default_backend()) where {T, TI}

    cell_mat = SMat{T}(cell)
    pbc_vec = SVec{Bool}(pbc)
    cutoff_T = T(cutoff)

    return _build_sorted_celllist(X, cell_mat, pbc_vec, cutoff_T, TI, backend)
end

"""
Internal implementation of sort-based cell list construction.
"""
function _build_sorted_celllist(X::AbstractVector{SVec{T}}, cell::SMat{T},
                                 pbc::SVec{Bool}, cutoff::T, ::Type{TI},
                                 backend) where {T, TI}
    nat = length(X)
    inv_cell, ncells_vec, lens = analyze_cell(cell, cutoff, TI)
    ncells_total = prod(ncells_vec)

    # Check for overflow
    if prod(BigInt.(ncells_vec)) > typemax(TI)
        error("""Ratio of simulation cell size to cutoff is very large.
                 Use a larger integer type (e.g. Int64), larger cutoff,
                 or smaller simulation cell.""")
    end

    # Step 1: Compute cell ID for each atom
    cell_ids = _compute_cell_ids(X, inv_cell, ncells_vec, pbc, TI)

    # Step 2: Sort atoms by cell ID
    perm = _get_sortperm(cell_ids, backend)
    sorted_cell_ids = cell_ids[perm]
    sorted_X = X[perm]

    # Step 3: Compute cell offsets (CSR-style)
    cell_offsets = _compute_cell_offsets(sorted_cell_ids, ncells_total, TI)

    return SortedCellList{T, TI, typeof(X), typeof(perm)}(
        sorted_X, X, perm, sorted_cell_ids, cell_offsets,
        cell, inv_cell, pbc, cutoff, ncells_vec, ncells_total
    )
end

"""
Compute cell ID for each atom position.
"""
function _compute_cell_ids(X::AbstractVector{SVec{T}}, inv_cell::SMat{T},
                           ncells::SVec{TI}, pbc::SVec{Bool}, ::Type{TI}) where {T, TI}
    nat = length(X)
    cell_ids = Vector{TI}(undef, nat)
    ns = ncells.data

    for i in 1:nat
        # Get cartesian cell index
        c = position_to_cell_index(inv_cell, X[i], ncells)
        # Apply boundary conditions
        c = bin_wrap_or_trunc.(c, pbc, ncells)
        # Convert to linear index
        cell_ids[i] = _sub2ind(ns, c)
    end

    return cell_ids
end

# Default sortperm for CPU
function _get_sortperm(cell_ids::AbstractVector{TI}, backend) where TI
    return TI.(sortperm(cell_ids))
end

# GPU extension point - override in NeighbourListsCUDAExt
function gpu_sortperm end

"""
Compute CSR-style offsets from sorted cell IDs.
`cell_offsets[c]:cell_offsets[c+1]-1` gives indices of atoms in cell `c`.
"""
function _compute_cell_offsets(sorted_cell_ids::AbstractVector{TI},
                               ncells_total::Integer, ::Type{TI}) where TI
    nat = length(sorted_cell_ids)
    # +1 for sentinel at end
    cell_offsets = zeros(TI, ncells_total + 1)

    if nat == 0
        cell_offsets .= one(TI)
        return cell_offsets
    end

    # Count atoms per cell
    for cid in sorted_cell_ids
        cell_offsets[cid + 1] += one(TI)
    end

    # Cumulative sum to get offsets (1-indexed)
    cell_offsets[1] = one(TI)
    for c in 1:ncells_total
        cell_offsets[c + 1] += cell_offsets[c]
    end

    return cell_offsets
end


# ==================== Lazy Neighbour Access for SortedCellList ====================

"""
    nsites(clist::SortedCellList)

Return the number of atoms/sites in the cell list.
"""
nsites(clist::SortedCellList) = length(clist.X)

"""
    cutoff(clist::SortedCellList)

Return the cutoff distance.
"""
cutoff(clist::SortedCellList) = clist.cutoff

"""
Get the range of adjacent cell indices to check for a given cell.
Calls `f(cell_index, shift_vector)` for each adjacent cell.

This version iterates directly without allocation for better performance.
For fully periodic systems, this is straightforward.
For non-periodic, we skip cells outside the domain bounds.
"""
@inline function _for_each_adjacent_cell(f::F, ci::SVec{TI}, ncells::SVec{TI},
                                          pbc::SVec{Bool}, nxyz::SVec{TI}) where {F, TI}
    ns = ncells.data

    # Compute the range of cells to visit in each dimension
    # For periodic: full range, wrapping handled by shift
    # For non-periodic: clamp to [1, ncells]
    x_lo = pbc[1] ? (ci[1] - nxyz[1]) : max(one(TI), ci[1] - nxyz[1])
    x_hi = pbc[1] ? (ci[1] + nxyz[1]) : min(ncells[1], ci[1] + nxyz[1])
    y_lo = pbc[2] ? (ci[2] - nxyz[2]) : max(one(TI), ci[2] - nxyz[2])
    y_hi = pbc[2] ? (ci[2] + nxyz[2]) : min(ncells[2], ci[2] + nxyz[2])
    z_lo = pbc[3] ? (ci[3] - nxyz[3]) : max(one(TI), ci[3] - nxyz[3])
    z_hi = pbc[3] ? (ci[3] + nxyz[3]) : min(ncells[3], ci[3] + nxyz[3])

    for cj_z in z_lo:z_hi
        cj_z_wrapped = pbc[3] ? bin_wrap(cj_z, ncells[3]) : cj_z
        shift_z = _compute_shift(cj_z, ncells[3], pbc[3])

        for cj_y in y_lo:y_hi
            cj_y_wrapped = pbc[2] ? bin_wrap(cj_y, ncells[2]) : cj_y
            shift_y = _compute_shift(cj_y, ncells[2], pbc[2])

            for cj_x in x_lo:x_hi
                cj_x_wrapped = pbc[1] ? bin_wrap(cj_x, ncells[1]) : cj_x
                shift_x = _compute_shift(cj_x, ncells[1], pbc[1])

                cj_wrapped = SVec{TI}(cj_x_wrapped, cj_y_wrapped, cj_z_wrapped)
                cell_shift = SVec{TI}(shift_x, shift_y, shift_z)

                linear_idx = _sub2ind(ns, cj_wrapped)
                f(linear_idx, cell_shift)
            end
        end
    end
end

# Legacy version that returns array (for compatibility) - avoid using in hot paths
function _adjacent_cells(ci::SVec{TI}, ncells::SVec{TI}, pbc::SVec{Bool},
                         nxyz::SVec{TI}) where TI
    result = Tuple{TI, SVec{TI}}[]
    seen = Set{Tuple{TI, SVec{TI}}}()

    _for_each_adjacent_cell(ci, ncells, pbc, nxyz) do linear_idx, cell_shift
        key = (linear_idx, cell_shift)
        if !(key in seen)
            push!(seen, key)
            push!(result, key)
        end
    end
    return result
end

@inline function _wrap_cell_index(i::TI, n::TI, pbc::Bool) where TI
    if pbc
        return bin_wrap(i, n)
    else
        return bin_trunc(i, n)
    end
end

@inline function _compute_shift(i::TI, n::TI, pbc::Bool) where TI
    if !pbc
        return zero(TI)
    elseif i <= 0
        return TI(-1)
    elseif i > n
        return TI(1)
    else
        return zero(TI)
    end
end

"""
    for_each_neighbour(f, clist::SortedCellList, i)

Call `f(j, R, shift)` for each neighbour j of atom i within cutoff.
`j` is the original (unsorted) atom index.
`R` is the displacement vector from i to j.
`shift` is the periodic cell shift.

This is designed to be efficient inside GPU kernels.
"""
@inline function for_each_neighbour(f::F, clist::SortedCellList{T,TI},
                                    i::Integer) where {F, T, TI}
    X_orig = clist.X_orig
    xi = X_orig[i]
    cutoff_sq = clist.cutoff^2
    inv_cell = clist.inv_cell
    ncells = clist.ncells
    pbc = clist.pbc
    cell_mat = clist.cell
    lens = lengths(cell_mat)
    cell_offsets = clist.cell_offsets
    perm = clist.perm

    # Get cell index of atom i
    ci = position_to_cell_index(inv_cell, xi, ncells)
    ci = bin_wrap_or_trunc.(ci, pbc, ncells)

    # How many cells to check in each direction
    nxyz = ceil.(TI, clist.cutoff * (ncells ./ abs.(lens)))

    # Iterate over adjacent cells (non-allocating version)
    _for_each_adjacent_cell(ci, ncells, pbc, nxyz) do cj_linear, cell_shift
        # Get atoms in this cell
        idx_start = cell_offsets[cj_linear]
        idx_end = cell_offsets[cj_linear + 1] - one(TI)

        for idx in idx_start:idx_end
            j_orig = perm[idx]  # Original atom index

            # Skip self-interaction (unless periodic image)
            if i == j_orig && all(cell_shift .== 0)
                continue
            end

            xj = X_orig[j_orig]

            # Compute displacement with periodic shift
            R = xj - xi + cell_mat' * cell_shift
            r_sq = dot(R, R)

            if r_sq < cutoff_sq
                f(j_orig, R, cell_shift)
            end
        end
    end
end

"""
    count_neighbours(clist::SortedCellList, i) -> Int

Count the number of neighbours of atom i within cutoff.
"""
function count_neighbours(clist::SortedCellList, i::Integer)
    count = 0
    for_each_neighbour(clist, i) do j, R, S
        count += 1
    end
    return count
end

"""
    get_neighbours(clist::SortedCellList, i) -> (js, Rs, Ss)

Get all neighbours of atom i. Returns vectors of indices, displacements, and shifts.
"""
function get_neighbours(clist::SortedCellList{T,TI}, i::Integer) where {T, TI}
    js = TI[]
    Rs = SVec{T}[]
    Ss = SVec{TI}[]

    for_each_neighbour(clist, i) do j, R, S
        push!(js, j)
        push!(Rs, R)
        push!(Ss, S)
    end

    return js, Rs, Ss
end


# ==================== Materialize PairList from SortedCellList ====================

"""
    materialize_pairlist(clist::SortedCellList; backend=CPU())

Convert a SortedCellList to a PairList by computing all pairs.

Uses a two-pass algorithm for parallelization:
1. Count neighbours per atom (parallel with KernelAbstractions)
2. Fill pair arrays (parallel with KernelAbstractions)

This approach enables multi-threaded CPU execution via KernelAbstractions.
"""
function materialize_pairlist(clist::SortedCellList{T, TI, <:Vector};
                              backend = default_backend()) where {T, TI}
    nat = nsites(clist)

    if nat == 0
        return PairList{T, TI, typeof(clist.X_orig), Vector{TI}, Vector{SVec{TI}}}(
            clist.X_orig, clist.cell, clist.cutoff,
            Vector{TI}(), Vector{TI}(), Vector{SVec{TI}}(), ones(TI, 1))
    end

    # Compute how many cells to check in each direction
    lens = lengths(clist.cell)
    nxyz = ceil.(TI, clist.cutoff * (clist.ncells ./ abs.(lens)))
    cutoff_sq = clist.cutoff^2

    # ===== Pass 1: Count neighbours per atom (parallel) =====
    counts = zeros(TI, nat)

    kernel1 = count_neighbours_kernel!(backend)
    kernel1(counts, clist.X_orig, clist.cell_offsets, clist.perm,
            clist.inv_cell, clist.cell, clist.ncells, clist.pbc,
            cutoff_sq, nxyz; ndrange=nat)
    synchronize(backend)

    # ===== Step 2: Compute offsets via prefix sum =====
    offsets = Vector{TI}(undef, nat + 1)
    offsets[1] = one(TI)
    for i in 1:nat
        offsets[i + 1] = offsets[i] + counts[i]
    end
    total_pairs = offsets[end] - one(TI)

    if total_pairs == 0
        first_arr = ones(TI, nat + 1)
        return PairList{T, TI, typeof(clist.X_orig), Vector{TI}, Vector{SVec{TI}}}(
            clist.X_orig, clist.cell, clist.cutoff,
            Vector{TI}(), Vector{TI}(), Vector{SVec{TI}}(), first_arr)
    end

    # ===== Step 3: Allocate pair arrays =====
    i_arr = Vector{TI}(undef, total_pairs)
    j_arr = Vector{TI}(undef, total_pairs)
    S_arr = Vector{SVec{TI}}(undef, total_pairs)

    # ===== Pass 2: Fill pair arrays (parallel) =====
    kernel2 = fill_pairs_kernel!(backend)
    kernel2(i_arr, j_arr, S_arr, offsets, clist.X_orig, clist.cell_offsets,
            clist.perm, clist.inv_cell, clist.cell, clist.ncells, clist.pbc,
            cutoff_sq, nxyz; ndrange=nat)
    synchronize(backend)

    return PairList{T, TI, typeof(clist.X_orig), Vector{TI}, Vector{SVec{TI}}}(
                    clist.X_orig, clist.cell, clist.cutoff,
                    i_arr, j_arr, S_arr, offsets)
end


# ==================== New API: PairList with backend ====================

"""
    PairList(X, cutoff, cell, pbc; backend=CPU(), int_type=Int32)

Create a PairList using the specified backend (CPU or GPU).

Uses sort-based cell list construction which is parallelizable.
"""
function PairList(X::AbstractVector{<:SVec}, cutoff::Real,
                  cell::AbstractMatrix, pbc;
                  backend = default_backend(),
                  int_type::Type = Int32)
    clist = build_cell_list(X, cutoff, cell, pbc;
                            int_type=int_type, backend=backend)
    return materialize_pairlist(clist; backend=backend)
end
