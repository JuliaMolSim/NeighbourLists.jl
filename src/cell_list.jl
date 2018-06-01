using Base.Threads

export npairs, nsites

PairList{T}(X::Vector{SVec{T}}, cutoff::AbstractFloat, cell::AbstractMatrix, pbc;
            int_type::Type = Int, store_first = true, sorted = true, fixcell = true) =
   _pairlist_(X, SMat{T}(cell), SVec{Bool}(pbc), T(cutoff), zero(int_type),
              store_first, sorted, fixcell)

PairList{T}(X::Matrix{T}, args...; kwargs...) =
   PairList(reinterpret(SVec{T}, X, (size(X,2),)), args...; varargs...)

npairs(nlist::PairList) = length(nlist.i)
nsites(nlist::PairList) = length(nlist.first) - 1


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
@inline function bin_trunc(i::Integer, n::Integer)
   if i <= 0;     i = 1
   elseif i > n;  i = n
   end
   return i
end

"apply bin_trunc only if open bdry"
@inline bin_trunc(i::Integer, pbc::Bool, n::Integer) = pbc ? i : bin_trunc(i, n)

"applies bin_trunc to open bdry and bin_wrap to periodic bdry"
@inline bin_wrap_or_trunc(i::Integer, pbc::Integer, n::Integer) =
      pbc ? bin_wrap(i, n) : bin_trunc(i, n)

"Map particle position to a (cartesian) cell index"
@inline position_to_cell_index{T, TI <: Integer}(
                        inv_cell::SMat{T}, x::SVec{T}, ns::SVec{TI}) =
   floor.(TI, ((inv_cell * x) .* ns + 1))


# ------------ The next two functions are the only dimension-dependent
#              parts of the code!

# an extension of sub2ind for the case when i is a vector (cartesian index)
@inline Base.sub2ind{TI <: Integer}(dims::NTuple{3,TI}, i::SVec{TI}) =
   sub2ind(dims, i[1], i[2], i[3])

lengths{T}(C::SMat{T}) =
   det(C) ./ SVec{T}(norm(C[2,:]×C[3,:]), norm(C[3,:]×C[1,:]), norm(C[1,:]×C[2,:]))


# --------------------------------------------------------------------------

function analyze_cell(cell, cutoff, _::TI) where TI
   # check the cell volume (allow only 3D volumes!)
   volume = abs(det(cell))
   @assert volume > 1e-12
   # precompute inverse of cell matrix for coordiate transformation
   inv_cell = inv(cell)
   # Compute distance of cell faces
   lens = lengths(cell)
   # Number of cells for cell subdivision
   ns_vec = max.(floor.(TI, lens / cutoff), 1)
   return inv_cell, ns_vec, lens
end

# multi-threading setup

function setup_mt(niter::TI, maxnt = MAX_THREADS[1]) where TI
   nt = minimum([6, nthreads(), ceil(TI, niter / 20), maxnt])
   nn = ceil.(TI, linspace(1, niter+1, nt+1))
   return nt, nn
end

# ==================== CellList Core  ================


function _celllist_(X::Vector{SVec{T}}, cell::SMat{T}, pbc::SVec{Bool},
            cutoff::T, _i::TI) where {T <: AbstractFloat, TI <: Integer}

   # ----- analyze cell -----
   nat = length(X)
   inv_cell, ns_vec, lens = analyze_cell(cell, cutoff, _i)
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
   last = Vector{TI}(ncells)
   next = Vector{TI}(nat)
   nats = zeros(TI, ncells)

   for i = 1:nat
      # Get cell index
      c = position_to_cell_index(inv_cell, X[i], ns_vec)
      # Periodic/non-periodic boundary conditions
      c = bin_wrap_or_trunc.(c, pbc, ns_vec)
      # linear cell index  # (+1 due to 1-based indexing)
      ci = sub2ind(ns, c)   # <<<<

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
   inv_cell, ns_vec, lens = analyze_cell(cell, cutoff, one(TI))

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
   absdist_t = Vector{T}[ Vector{T}()  for n = 1:nt ]  # r_ij ~ norm(X[i]-X[j])
   distvec_t = Vector{SVec{T}}[ Vector{SVec{T}}()  for n = 1:nt ]  # ~ X[i] - X[j]
   # give size hints
   sz = (nat ÷ nt + nt) * nneigs_guess
   for n = 1:nt
      sizehint!(first_t[n], sz)
      sizehint!(secnd_t[n], sz)
      sizehint!(absdist_t[n], sz)
      sizehint!(distvec_t[n], sz)
   end

   # We need the shape of the bin ( bins[:, i] = cell[i,:] / ns[i] )
   bins = cell' ./ ns_vec

   # Find out over how many neighbor cells we need to loop (if the box is small)
   nxyz = ceil.(TI, cutoff * (ns_vec ./ lens))
   cxyz = CartesianIndex(nxyz.data)
   xyz_range = CartesianRange(- cxyz, cxyz)

   # Loop over threads
   @threads for it = 1:nt
      for i = nn[it]:(nn[it+1]-1)
         # current atom position
         _find_neighbours_!(i, clist, ns_vec, bins, xyz_range,
                     first_t[it], secnd_t[it], absdist_t[it], distvec_t[it])
      end # for i = 1:nat
   end # @threads

   # piece them back together
   sz = sum( length(first_t[i]) for i = 1:nt )
   first = first_t[1];     sizehint!(first, sz)
   secnd = secnd_t[1];     sizehint!(secnd, sz)
   absdist = absdist_t[1]; sizehint!(absdist, sz)
   distvec = distvec_t[1]; sizehint!(distvec, sz)
   for it = 2:nt
      append!(first, first_t[it])
      append!(secnd, secnd_t[it])
      append!(absdist, absdist_t[it])
      append!(distvec, distvec_t[it])
   end

   # Build return tuple
   return first, secnd, absdist, distvec
end



function _find_neighbours_!(i, clist, ns_vec::SVec{TI}, bins, xyz_range,
                            first, secnd, absdist, distvec) where TI
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
   dxi = xi - bins * (ci - 1)

   # Apply periodic boundary conditions as well now
   ci = bin_wrap_or_trunc.(ci0, pbc, ns_vec)

   for ixyz in xyz_range
      # convert cartesian index to SVector
      xyz = SVec{TI}(ixyz.I)
      # get the bin index
      cj = bin_wrap.(ci + xyz, pbc, ns_vec)
      # skip this bin if not inside the domain
      all(1 .<= cj .<= ns_vec) || continue
      # linear cell index
      ncj = sub2ind(ns_vec.data, cj)
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
            dxj = xj - bins * (cj - 1)
            # Compute distance between atoms
            dx = dxj - dxi + off
            norm_dx_sq = dx ⋅ dx

            # append to the list
            if norm_dx_sq < cutoff_sq
               push!(first, i)
               push!(secnd, j)
               push!(distvec, dx)
               push!(absdist, sqrt(norm_dx_sq))
            end
         end  # if i != j || any(xyz .!= 0)

         # go to the next atom in the current cell
         j = next[j];
      end # while j > 0 (loop over atoms in current cell)
   end # loop over neighbouring bins
end


function _pairlist_(X::Vector{SVec{T}}, cell::SMat{T}, pbc::SVec{Bool},
            cutoff::T, _i::TI, store_first::Bool, sorted::Bool, fixcell::Bool) where {T, TI}

   # temporary (?) fix to make sure all atoms are within the cell
   if fixcell
      X, cell = _fix_cell_(X, cell, pbc)
   end

   clist = _celllist_(X, cell, pbc, cutoff, _i)
   i, j, r, R = _pairlist_(clist)

   if store_first
      first = get_first(i, length(X))
   else
      first = zeros(int_type, 0)
   end

   if store_first && sorted
      sort_neigs!(j, r, R, first)
   end

   return PairList(X, cutoff, i, j, r, R, first)
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
function get_first{TI}(i::Vector{TI}, nat::Integer = i[end])
   # compute the first index for each site
   first = Vector{TI}(nat + 1)
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
   first[n:end] = length(i)+1
   return first
end


"""
`sort_neigs!(j, r, R, first)`

sorts each sub-range of `j` corresponding to one site  in ascending order
and applies the same permutation to `r, R, S`.
"""
function sort_neigs!(j, r, R, first)
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
               r[rg] = r[rg_perm]
               R[rg] = R[rg_perm]
            end
         end
      end
   end
end

"""
`max_neigs(nlist::PairList) -> Integer`

returns the maximum number of neighbours that any atom in the
neighbourlist has.
"""
function max_neigs(nlist::PairList)
      maximum( nlist.first[n+1]-nlist.first[n]
               for n = 1:(length(nlist.first)-1) )
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
      # we need to multiply the cell by the correct lambda
      C = hcat( max_lam[1] * C[1,:], max_lam[2] * C[2,:], max_lam[3] * C[3,:])'
   end
   return X, C
end
