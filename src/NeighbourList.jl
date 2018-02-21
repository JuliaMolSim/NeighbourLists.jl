module NeighbourList

using StaticArrays

const SVec{T} = SVector{3, T}
const SMat{T} = SMatrix{3, 3, T, 9}


# ====================== Cell Index Algebra =====================



"""
Map i back to the interval [0,n) by shifting by integer multiples of n
"""
@inline function bin_wrap(i::Integer, n::Integer)
   while i <= 0;  i += n; end
   while i > n;  i -= n; end
   return i;
end


"""
Map i back to the interval [0,n) by assigning edge value if outside interval
"""
@inline function bin_trunc(i::Integer, n::Integer)
   if i <= 0
      i = 1
   elseif i > n
      i = n
   end
   return i
end


@inline bin_trunc{TI}(cj::SVec{TI}, pbc::SVec{Bool}, ns::SVec{TI}) =
   SVec{TI}(pbc[1] ? cj[1] : bin_trunc(cj[1], ns[1]),
            pbc[2] ? cj[2] : bin_trunc(cj[2], ns[2]),
            pbc[3] ? cj[3] : bin_trunc(cj[3], ns[3]) )

@inline bin_wrap_or_trunc{TI}(cj::SVec{TI}, pbc::SVec{Bool}, ns::SVec{TI}) =
   SVec{TI}(pbc[1] ? bin_wrap(cj[1], ns[1]) : bin_trunc(cj[1], ns[1]),
            pbc[2] ? bin_wrap(cj[2], ns[2]) : bin_trunc(cj[2], ns[2]),
            pbc[3] ? bin_wrap(cj[3], ns[3]) : bin_trunc(cj[3], ns[3]) )


"""
Map particle position to a cell index
"""
@inline function position_to_cell_index{T, TI <: Integer}(
                        inv_cell::SMat{T}, x::SVec{T}, ns::SVec{TI})
   y = inv_cell * x
   return SVec{TI}( floor(TI, y[1]*ns[1]+1),
                    floor(TI, y[2]*ns[2]+1),
                    floor(TI, y[3]*ns[3]+1) )
end


# ==================== Actual Neighbour List Construction ================

# TODO: consider redefining the cell in terms of the
#       extent of X

function neighbour_list{T}(cell::SMat{T}, pbc::SVec{Bool}, X::Vector{SVec{T}},
                           cutoff::T, int_type::Type = Int)
   # @assert int_type <: Integer
   return _neighbour_list_(cell, pbc, X, cutoff, zero(int_type))
end


"""
TODO: write documentation
"""
function _neighbour_list_{T, TI}(cell::SMat{T}, pbc::SVec{Bool}, X::Vector{SVec{T}},
                           cutoff::T, ::TI)

   # ----------- analyze cell --------------
   # check the cell volume (allow only 3D volumes!)
   volume = det(cell)
   if volume < 1e-12
      error("Zero cell volume.")
   end
   # some cell geometry
   inv_cell = inv(cell)
   norm1 = cell[2, :] × cell[3, :]    # face normals
   norm2 = cell[3, :] × cell[1, :]
   norm3 = cell[1, :] × cell[2, :]
   # Compute distance of cell faces
   len1 = volume / norm(norm1);
   len2 = volume / norm(norm2);
   len3 = volume / norm(norm3);
   # Number of cells for cell subdivision
   n1 = max(floor(TI, len1 / cutoff), 1);
   n2 = max(floor(TI, len2 / cutoff), 1);
   n3 = max(floor(TI, len3 / cutoff), 1);
   ns = tuple(n1, n2, n3)
   ns_vec = SVec(n1, n2, n3)

   if BigInt(n1) * BigInt(n2) * BigInt(n3) > typemax(TI)
      error("""Ratio of simulation cell size to cutoff is very large.
               Are you using a cell with lots of vacuum? To fix this
               use a larger integer type (e.g. Int128) or a
               larger cut-off.""")
   end

   # just to double-check everything is sane?
   @assert all(ns .> 0)

   # Find out over how many neighbor cells we need to loop (if the box is small
   # TODO PBC
   nx = ceil(TI, cutoff * n1 / len1)
   ny = ceil(TI, cutoff * n2 / len2)
   nz = ceil(TI, cutoff * n3 / len3)

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
      c = bin_wrap_or_trunc(c, pbc, ns_vec)
      # linear cell index  # (+1 due to 1-based indexing)
      ci = sub2ind(ns, c[1], c[2], c[3])

      @assert all(1 .<= c .<= ns_vec)
      @assert 1 .<= ci .<= ncells

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
   first   = Vector{TI}();      sizehint!(first, szhint)  # i
   secnd   = Vector{TI}();      sizehint!(secnd, szhint)  # j -> (i,j) is a bond
   absdist = Vector{T}();       sizehint!(absdist, szhint) # Xj - Xi
   distvec = Vector{SVec{T}}(); sizehint!(distvec, szhint) # r_ij
   shift   = Vector{SVec{T}}(); sizehint!(shift, szhint)   # cell shifts

   # funnily testing with cutoff^2 actually makes a measurable difference
   # which suggests we are pretty close to the limit
   cutoff_sq = cutoff^2

   # We need the shape of the bin
   bin1 = cell[1, :] / n1
   bin2 = cell[2, :] / n2
   bin3 = cell[3, :] / n3

   # Loop over atoms
   for i = 1:nat
      # current atom position
      xi = X[i]
      # cell index (cartesian) of xi
      ci0 = position_to_cell_index(inv_cell, xi, ns_vec)

      # Truncate if non-periodic and outside of simulation domain
      # here, we don't yet want to wrap the pbc as well
      ci = bin_trunc(ci0, pbc, ns_vec)
      # dri is the position relative to the lower left corner of the bin
      dxi = xi - (ci[1]-1) * bin1 - (ci[2]-1) * bin2 - (ci[3]-1) * bin3

      # Apply periodic boundary conditions as well
      # (technically we only need to wrap here, since we've already truncated
      #  but it makes very little difference for performance)
      ci = bin_wrap_or_trunc(ci0, pbc, ns_vec)

      # Loop over neighbouring bins
      for z = -nz:nz
         cj3 = ci[3] + z
         if pbc[3]
            cj3 = bin_wrap(cj3, n3)
         end
         # Skip to next z value if neighbour cell is out of simulation bounds
         if cj3 <= 0 || cj3 > n3
            continue
         end

         for y = -ny:ny
            cj2 = ci[2] + y
            if pbc[2]
               cj2 = bin_wrap(cj2, n2)
            end
            # Skip to next y value if cell is out of simulation bounds
            if cj2 <= 0 || cj2 > n2
               continue
            end

            for x = -nx:nx
               # Bin index of neighbouring bin
               cj1 = ci[1] + x
               if pbc[1]
                  cj1 = bin_wrap(cj1, n1)
               end
               # Skip to next x value if cell is out of simulation bounds
               if cj1 <= 0 || cj1 > n1
                  continue
               end

               # linear cell index
               ncj = sub2ind(ns, cj1, cj2, cj3)

               # Offset of the neighboring bins
               off = x * bin1 + y * bin2 + z * bin3

               # Loop over all atoms in neighbouring bin (all potential
               # neighbours in the bin with linear index cj1)
               j = seed[ncj] # the first atom in the ncj cell
               while j > 0
                  if i != j || x != 0 || y != 0 || z != 0
                     xj = X[j] # position of current neighbour

                     # we need to find the cell index again, because this is
                     # not really the cell index, but it could be outside
                     # the domain -> i.e. this only makes a difference for pbc
                     # TODO: which suggests there is an optimisation to be done
                     cj = position_to_cell_index(inv_cell, xj, ns_vec)
                     cj = bin_trunc(cj, pbc, ns_vec)

                     # drj is position relative to lower left corner of the bin
                     dxj = xj - (cj[1]-1) * bin1 - (cj[2]-1) * bin2 - (cj[3]-1) * bin3

                     # Compute distance between atoms
                     dx = dxj - dxi + off
                     norm_dx_sq = dot(dx, dx)

                     if norm_dx_sq < cutoff_sq
                        push!(first, i)
                        push!(secnd, j)
                        push!(distvec, dx)
                        push!(absdist, sqrt(norm_dx_sq))
                        push!(shift, SVec{TI}((ci0[1] - cj[1] + x) ÷ n1,
                                              (ci0[2] - cj[2] + y) ÷ n2,
                                              (ci0[3] - cj[3] + z) ÷ n3))
                     end
                  end  # if i != j || x != 0 || y != 0 || z != 0

                  # go to the next atom in the current cell
                  j = next[j];

               end # while j > 0
            end # for x
         end # for y
      end # for z
   end # for i = 1:nat

   # Build return tuple
   return first, secnd, absdist, distvec, shift
end


end # module
