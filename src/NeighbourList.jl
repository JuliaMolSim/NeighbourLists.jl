module NeighbourList

using StaticArrays

const SVec{T} = SVector{3, T}
const SMat{T} = SMatrix{3, 3, T, 9}


# ====================== Cell Index Algebra =====================

# TODO: should indexing of the cells start at 0?


"""
Map i back to the interval [0,n) by shifting by integer multiples of n
"""
@inline function bin_wrap(i::Integer, n::Integer)
   while i <= 0;  i += n; end
   while i > n;  i -= n; end
   return i;
end
# TODO: redo this using round?

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


# matrix-vector multiplication
# mat_mul_vec(double *mat, double *vin, double *vout)

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

# position_to_cell_index(double *inv_cell, double *ri, int n1, int n2, int n3,
#                        int *c1, int *c2, int *c3)
#  >>>> The output was stored in c1, c2, c3!


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
   inv_cell = inv(cell)
   cell1 = cell[1, :]       # cell vectors
   cell2 = cell[2, :]
   cell3 = cell[3, :]
   norm1 = cell2 × cell3    # face normals
   norm2 = cell3 × cell1
   norm3 = cell1 × cell2
   # check the cell volume (allow only 3D volumes!)
   volume = dot(cell3, norm3)
   if volume < 1e-12
      error("Zero cell volume.")
   end
   # normalise the face normals
   len1, len2, len3 = norm(norm1), norm(norm2), norm(norm3)
   norm1 *= volume/len1
   norm2 *= volume/len2
   norm3 *= volume/len3
   # Compute distance of cell faces
   len1 = volume/len1;
   len2 = volume/len2;
   len3 = volume/len3;
   # Number of cells for cell subdivision
   n1 = max(floor(TI, len1 / cutoff), 1);
   n2 = max(floor(TI, len2 / cutoff), 1);
   n3 = max(floor(TI, len3 / cutoff), 1);
   ns = tuple(n1, n2, n3)
   ns_vec = SVec(n1, n2, n3)

   if BigInt(n1) * BigInt(n2) * BigInt(n3) > typemax(Int)
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
      # TODO PBC
      c1 = pbc[1] ? bin_wrap(c[1], n1) : bin_trunc(c[1], n1)
      c2 = pbc[2] ? bin_wrap(c[2], n2) : bin_trunc(c[2], n2)
      c3 = pbc[3] ? bin_wrap(c[3], n3) : bin_trunc(c[3], n3)

      # linear cell index  # (+1 due to 1-based indexing)
      ci = sub2ind(ns, c1, c2, c3)

      # TODO: rewrite more nicely :)   1 <= c1 = n1 etc
      #       or  all(1 .<= cs .<= ns)
      @assert c1 > 0 && c1 <= n1
      @assert c2 > 0 && c2 <= n2
      @assert c3 > 0 && c3 <= n3
      @assert ci > 0 && ci <= ncells

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
   # Neighbour list counter and size
   szhint = nat*6    # Initial guess for neighbour list size

   # neighbourlist information
   first   = Vector{TI}();      sizehint!(first, szhint)  # i
   secnd   = Vector{TI}();      sizehint!(secnd, szhint)  # j -> (i,j) is a bond
   absdist = Vector{T}();       sizehint!(absdist, szhint)       # Xj - Xi
   distvec = Vector{SVec{T}}(); sizehint!(distvec, szhint) # r_ij
   shift   = Vector{SVec{T}}(); sizehint!(shift, szhint)     # cell shifts

   # funnily testing with cutoff^2 actually makes a measurable difference
   cutoff_sq = cutoff^2

   # We need the shape of the bin
   bin1 = cell1 / n1
   bin2 = cell2 / n2
   bin3 = cell3 / n3

   # Loop over atoms
   for i = 1:nat
      xi = X[i]
      ci0 = position_to_cell_index(inv_cell, xi, ns_vec)

      # Truncate if non-periodic and outside of simulation domain
      ci1 = pbc[1]  ?  ci0[1]  :  bin_trunc(ci0[1], n1)
      ci2 = pbc[2]  ?  ci0[2]  :  bin_trunc(ci0[2], n2)
      ci3 = pbc[3]  ?  ci0[3]  :  bin_trunc(ci0[3], n3)

      # dri is the position relative to the lower left corner of the bin
      dxi = xi - (ci1-1) * bin1 - (ci2-1) * bin2 - (ci3-1) * bin3

      # Apply periodic boundary conditions
      # TODO PBC; ALSO why is bin_trunc run a second time???
      ci1 = pbc[1]  ?  bin_wrap(ci0[1], n1)  :  bin_trunc(ci0[1], n1)
      ci2 = pbc[2]  ?  bin_wrap(ci0[2], n2)  :  bin_trunc(ci0[2], n2)
      ci3 = pbc[3]  ?  bin_wrap(ci0[3], n3)  :  bin_trunc(ci0[3], n3)

      # Loop over neighbouring bins
      for z = -nz:nz
         cj3 = ci3 + z
         if pbc[3]    # TODO PBC
            cj3 = bin_wrap(cj3, n3)
         end

         # Skip to next z value if neighbour cell is out of simulation bounds
         if cj3 <= 0 || cj3 > n3
            continue
         end

         # cj3 = bin_trunc(cj3, n3)  #  THIS SHOULD NOT BE NECESSARY?
         ncj3 = n2 * (cj3-1)    # WHAT DOES THIS DO?
         off3 = z * bin3

         for y = -ny:ny
            cj2 = ci2 + y
            if pbc[2]   # TODO PBC
               cj2 = bin_wrap(cj2, n2)
            end

            # Skip to next y value if cell is out of simulation bounds
            if cj2 <= 0 || cj2 > n2
               continue
            end

            # cj2 = bin_trunc(cj2, n2)  SHOUDL NOT BE NECESSARY
            ncj2 = n1 * (cj2-1 + ncj3)   # what does this do?
            off2 = off3 + y * bin2;

            for x = -nx:nx
               # Bin index of neighbouring bin
               cj1 = ci1 + x
               if pbc[1]   # TODO PBC
                  cj1 = bin_wrap(cj1, n1)
               end

               # Skip to next x value if cell is out of simulation bounds
               if cj1 <= 0 || cj1 > n1
                  continue
               end

               # cj1 = bin_trunc(cj1, n1)   should not be necessary
               ncj = cj1 + ncj2   # linear cell-index of potential neighbour j
               @assert ncj == sub2ind(ns, cj1, cj2, cj3)
               # TODO: switch to sub2ind??

               # Offset of the neighboring bins
               off = off2 + x * bin1

               # Loop over all atoms in neighbouring bin (all potential
               # neighbours in the bin with linear index cj1)
               j = seed[ncj] # the first atom in the ncj cell
               while j > 0
                  if i != j || x != 0 || y != 0 || z != 0
                     xj = X[j] # position of current neighbour

                     # TODO: Why do we need to find the cell index again?
                     cj = position_to_cell_index(inv_cell, xj, ns_vec)
                     cj = MVector{3, TI}(cj_)

                     # Truncate if non-periodic and outside of simulation domain
                     if !pbc[1];  cj[1] = bin_trunc(cj[1], n1);  end
                     if !pbc[2];  cj[2] = bin_trunc(cj[2], n2);  end
                     if !pbc[3];  cj[3] = bin_trunc(cj[3], n3);  end

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
