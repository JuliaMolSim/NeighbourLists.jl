
# ================= types and wrapper stuff ======================

"""
`CellList` stores a neighbourlist that is constructed from a cell-list
(it isn't actually stored as a cell-list).
Typically, this is constructed using
```
CellList(X, cutoff, cell, pbc)
```
where

* `X` : positions, either as 3 x N matrix or as `Vector{SVec}`
* `cutoff` : positive real value
* `cell` : 3 x 3 matrix, with rows denoting the cell vectors
* `pbc` : 3-tuple or array, storing periodicity information

### Kw-args:

* `int_type` : default is `Int`
* `store_first` : whether to store the array of first indices, default `true`
* `sorted` : whether to sort the `j` vector, default `false`

### CellList fields

`i, j, r, R, first`, where

`(i[n], j[n])` denotes the indices of a neighbour pair, `r[n]` the distance
between those atoms, `R[n]` the vectorial distance, note this is identical to
`X[i[n]]-X[j[n]]` without periodic b.c.s, but with periodic boundary conditions
it is different. `first[m]` contains the index to the first `(i[n], j[n])` for
which `i[n] == first[m]`, i.e., `(j, first)` essentially defines a compressed
columns storage of the adjacancy matrix.
"""
struct CellList{T <: AbstractFloat, TI <: Integer}
   X::Vector{SVec{T}}
   cutoff::T
   i::Vector{TI}
   j::Vector{TI}
   r::Vector{T}
   R::Vector{SVec{T}}
   S::Vector{SVec{TI}}
   first::Vector{TI}
end

# TODO: consider redefining the cell in terms of the
#       extent of X

function CellList{T}(X::Vector{SVec{T}}, cutoff::AbstractFloat,
                     cell::AbstractMatrix, pbc;
                     int_type::Type = Int, store_first = true,
                     sorted = false, neig_guess = 12)
   i, j, r, R, S = _neighbour_list_(SMat{T}(cell...), SVec{Bool}(pbc...), X,
                                    cutoff, zero(int_type), neig_guess)
   if store_first
      first = get_first(i, length(X))
   else
      first = zeros(int_type, 0)
   end

   if store_first && sorted
      sort_neigs!(j, r, R, S, first)
   end

   return CellList(X, cutoff, i, j, r, R, S, first)
end

CellList{T}(X::Matrix{T}, args...; kwargs...) =
      CellList(reinterpret(SVec{T}, X, (size(X,2),)), args...; varargs...)

npairs(nlist::CellList) = length(nlist.i)
nsites(nlist::CellList) = length(nlist.first) - 1


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



# ==================== CellList Core  ================


"""
_neighbour_list_(cell, pbc, X, cutoff, _): inner neighbourlist assembly

*   `cell::SMatrix{D, D, T}`: rows are the cell vectors
*    `pbc::SVector{3, Bool}`: flags for periodic bcs
*      `X::Vector{SVec{T}}` : positions
* `cutoff::T`               : cutoff radius
* `     _::TI`              : a number specifying the integer type to be used
"""
function _neighbour_list_{T, TI}(cell::SMat{T}, pbc::SVec{Bool}, X::Vector{SVec{T}},
                           cutoff::T, ::TI, neigs_guess::Integer = 12)

   # ----------- analyze cell --------------
   # check the cell volume (allow only 3D volumes!)
   volume = det(cell)
   if volume < 1e-12
      error("(near) Zero cell volume.")
   end
   # precompute inverse of cell matrix for coordiate transformation
   inv_cell = inv(cell)
   # Compute distance of cell faces
   lens = lengths(cell)
   # Number of cells for cell subdivision
   ns_vec = max.(floor.(TI, lens / cutoff), 1)
   ns = ns_vec.data   # a tuple

   if prod(BigInt.(ns_vec)) > typemax(TI)
      error("""Ratio of simulation cell size to cutoff is very large.
               Are you using a cell with lots of vacuum? To fix this
               use a larger integer type (e.g. Int128), a
               larger cut-off, or a smaller simulation cell.""")
   end

   # Find out over how many neighbor cells we need to loop (if the box is small)
   nxyz = ceil.(TI, cutoff * (ns_vec ./ lens))
   cxyz = CartesianIndex(nxyz.data)
   xyz_range = CartesianRange(- cxyz, cxyz)

   # ------------ Sort particles into bins -----------------
   nat = length(X)
   ncells = prod(ns_vec)
   # data structure to store a linked list for each bin
   seed = fill(TI(-1), ncells)
   last = Vector{TI}(ncells)
   next = Vector{TI}(nat)

   for i = 1:nat
      # Get cell index
      c = position_to_cell_index(inv_cell, X[i], ns_vec)
      # Periodic/non-periodic boundary conditions
      c = bin_wrap_or_trunc.(c, pbc, ns_vec)
      # linear cell index  # (+1 due to 1-based indexing)
      ci = sub2ind(ns, c)   # <<<<
      # sanity check
      # @assert all(1 .<= c .<= ns_vec)
      # @assert 1 <= ci <= ncells

      # Put atom into appropriate bin (linked list)
      if seed[ci] < 0
         next[i] = -1;
         seed[ci] = i;
         last[ci] = i;
      else
         next[i] = -1;
         next[last[ci]] = i;
         last[ci] = i;
      end
   end

   # ------------ Start actual neighbourlist assembly ----------
   # allocate neighbourlist information (can make a better guess?)
   szhint = nat*12    # Initial guess for neighbour list size
   first   = Vector{TI}();       sizehint!(first, szhint)   # i
   secnd   = Vector{TI}();       sizehint!(secnd, szhint)   # j -> (i,j) is a bond
   absdist = Vector{T}();        sizehint!(absdist, szhint) # Xj - Xi
   distvec = Vector{SVec{T}}();  sizehint!(distvec, szhint) # r_ij
   shift   = Vector{SVec{TI}}(); sizehint!(shift, szhint)   # cell shifts

   # funnily testing with cutoff^2 actually makes a measurable difference
   # which suggests we are pretty close to the performance limit
   cutoff_sq = cutoff^2

   # We need the shape of the bin ( bins[:, i] = cell[i,:] / ns[i] )
   bins = cell' ./ ns_vec

   # Loop over atoms
   for i = 1:nat
      # current atom position
      xi = X[i]
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
         ncj = sub2ind(ns, cj)    # <<<<<<<<
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
                  # push!(shift, (ci0 - cj + xyz) .÷ ns_vec)
               end
            end  # if i != j || any(xyz .!= 0)

            # go to the next atom in the current cell
            j = next[j];
         end # while j > 0 (loop over atoms in current cell)
      end # loop over neighbouring bins
   end # for i = 1:nat

   # Build return tuple
   return first, secnd, absdist, distvec, shift
end



# ==================== Bonus Stuff =========================

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
function sort_neigs!(j, r, R, S, first)
   for n = 1:length(first)-1
      if first[n+1] > first[n] + 1
         rg = first[n]:first[n+1]-1
         I = sortperm(@view j[rg])
         rg_perm = rg[I]
         j[rg] = j[rg_perm]
         r[rg] = r[rg_perm]
         R[rg] = R[rg_perm]
         # S[rg] = S[rg_perm]
      end
   end
end
