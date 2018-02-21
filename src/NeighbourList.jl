module NeighbourList

using StaticArrays, ResumableFunctions

const SVec{T} = SVector{3, T}
const SMat{T} = SMatrix{3, 3, T, 9}


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

# an extension of sub2ind for the case when i is a vector (cartesian index)
@inline Base.sub2ind{TI <: Integer}(dims::NTuple{3,TI}, i::SVec{TI}) =
   sub2ind(dims, i[1], i[2], i[3])


lengths{T}(C::SMat{T}) =
   det(C) ./ SVec{T}(norm(C[2,:]×C[3,:]), norm(C[3,:]×C[1,:]), norm(C[1,:]×C[2,:]))


# ==================== GateWay Routines  ================

# TODO: consider redefining the cell in terms of the
#       extent of X

function neighbour_list{T}(cell::SMat{T}, pbc::SVec{Bool}, X::Vector{SVec{T}},
                           cutoff::T, int_type::Type = Int)
   return _neighbour_list_(cell, pbc, X, cutoff, zero(int_type))
end


# ==================== Main NeighbourList Core  ================


"""
_neighbour_list_(cell, pbc, X, cutoff, _): inner neighbourlist assembly

*   `cell::SMatrix{D, D, T}`: rows are the cell vectors
*    `pbc::SVector{3, Bool}`: flags for periodic bcs
*      `X::Vector{SVec{T}}` : positions
* `cutoff::T`               : cutoff radius
* `     _::TI`              : a number specifying the integer type to be used
"""
function _neighbour_list_{T, TI}(cell::SMat{T}, pbc::SVec{Bool}, X::Vector{SVec{T}},
                           cutoff::T, _::TI)

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

   # Find out over how many neighbor cells we need to loop (if the box is small
   nxyz = ceil.(TI, cutoff * (ns_vec ./ lens))
   cxyz = CartesianIndex(nxyz.data)
   xyz_range = CartesianRange(- cxyz, cxyz)

   # ------------ Sort particles into bins -----------------
   nat = length(X)
   ncells = prod(ns)
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
      ci = sub2ind(ns, c)
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
   szhint = nat*6    # Initial guess for neighbour list size
   first   = Vector{TI}();      sizehint!(first, szhint)   # i
   secnd   = Vector{TI}();      sizehint!(secnd, szhint)   # j -> (i,j) is a bond
   absdist = Vector{T}();       sizehint!(absdist, szhint) # Xj - Xi
   distvec = Vector{SVec{T}}(); sizehint!(distvec, szhint) # r_ij
   shift   = Vector{SVec{T}}(); sizehint!(shift, szhint)   # cell shifts

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
         ncj = sub2ind(ns, cj)
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
                  push!(shift, (ci0 - cj + xyz) .÷ ns_vec)
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


end # module
