using NearestNeighbors, Distances

import Distances: evaluate

export nn_list

#
# function NNList{T}(X::Vector{SVec{T}}, cutoff:T,
#                         cell::SMat{T}, pbc::SVec{Bool}; int_type = Int)
#
# end


immutable Periodic{T} <: Metric
   L::SVec{T}
   pbc::SVec{Bool}
end

function evaluate{T}(d::Periodic{T}, x::AbstractVector, y::AbstractVector)
   s = 0.0
   for i = 1:3
      if d.pbc[i]
         m = mod(x[i]-y[i], d.L[i])
         s += min( m, d.L[i]-m )^2
      else
         s += (x[i]-y[i])^2
      end
   end
   return sqrt(s)
end

function nn_list{T, TI}(X::Vector{SVec{T}}, cutoff::T,
                        cell::SMat{T}, pbc::SVec{Bool}, _::TI)
   @assert cell == SMat{T}(diagm(diag(cell)))   # require cubic cell

   # construct a periodic metric
   d = Periodic(SVec{T}(diag(cell)), pbc)

   # construct the ball tree
   tree = BallTree(X, d)

   i = TI[]
   j = TI[]
   r = T[]
   R = SVec{T}[]

   for n = 1:length(X)
      neigs = inrange(tree, X[n], cutoff, true)
      for a = 1:length(neigs)
         if neigs[a] != n
            push!(i, n)
            push!(j, neigs[a])
         end 
      end
   end

   return i, j, r, R
end




# function _nn_list_{T, TI}(cell::SMat{T}, pbc::SVec{Bool}, X::Vector{SVec{T}},
#                            cutoff::T, _::TI)
#
#    # ----------- analyze cell --------------
#    # check the cell volume (allow only 3D volumes!)
#    volume = det(cell)
#    if volume < 1e-12
#       error("(near) Zero cell volume.")
#    end
#    # precompute inverse of cell matrix for coordiate transformation
#    inv_cell = inv(cell)
#    # Compute distance of cell faces
#    lens = lengths(cell)
#    # Number of cells for cell subdivision
#    ns_vec = max.(floor.(TI, lens / cutoff), 1)
#    ns = ns_vec.data   # a tuple
#
#    if prod(BigInt.(ns_vec)) > typemax(TI)
#       error("""Ratio of simulation cell size to cutoff is very large.
#                Are you using a cell with lots of vacuum? To fix this
#                use a larger integer type (e.g. Int128), a
#                larger cut-off, or a smaller simulation cell.""")
#    end
#
#    # Find out over how many neighbor cells we need to loop (if the box is small)
#    nxyz = ceil.(TI, cutoff * (ns_vec ./ lens))
#    cxyz = CartesianIndex(nxyz.data)
#    xyz_range = CartesianRange(- cxyz, cxyz)
#
#    # ------------ Create a NearestNeighbours.jl Tree -----------------
#    nat = length(X)
#
#
#    # ------------ Start actual neighbourlist assembly ----------
#    # allocate neighbourlist information (can make a better guess?)
#    szhint = nat*6    # Initial guess for neighbour list size
#    first   = Vector{TI}();       sizehint!(first, szhint)   # i
#    secnd   = Vector{TI}();       sizehint!(secnd, szhint)   # j -> (i,j) is a bond
#    absdist = Vector{T}();        sizehint!(absdist, szhint) # Xj - Xi
#    distvec = Vector{SVec{T}}();  sizehint!(distvec, szhint) # r_ij
#    shift   = Vector{SVec{TI}}(); sizehint!(shift, szhint)   # cell shifts
#
#    # funnily testing with cutoff^2 actually makes a measurable difference
#    # which suggests we are pretty close to the performance limit
#    cutoff_sq = cutoff^2
#
#    # We need the shape of the bin ( bins[:, i] = cell[i,:] / ns[i] )
#    bins = cell' ./ ns_vec
#
#    # Loop over atoms
#    for i = 1:nat
#       # current atom position
#       xi = X[i]
#       # cell index (cartesian) of xi
#       ci0 = position_to_cell_index(inv_cell, xi, ns_vec)
#
#       # Truncate if non-periodic and outside of simulation domain
#       # (here, we don't yet want to wrap the pbc as well)
#       ci = bin_trunc.(ci0, pbc, ns_vec)
#       # dxi is the position relative to the lower left corner of the bin
#       dxi = xi - bins * (ci - 1)
#
#       # Apply periodic boundary conditions as well now
#       ci = bin_wrap_or_trunc.(ci0, pbc, ns_vec)
#
#       for ixyz in xyz_range
#          # convert cartesian index to SVector
#          xyz = SVec{TI}(ixyz.I)
#          # get the bin index
#          cj = bin_wrap.(ci + xyz, pbc, ns_vec)
#          # skip this bin if not inside the domain
#          all(1 .<= cj .<= ns_vec) || continue
#          # linear cell index
#          ncj = sub2ind(ns, cj)    # <<<<<<<<
#          # Offset of the neighboring bins
#          off = bins * xyz
#
#          # Loop over all atoms in neighbouring bin (all potential
#          # neighbours in the bin with linear index cj1)
#          j = seed[ncj] # the first atom in the ncj cell
#          while j > 0
#             if i != j || any(xyz .!= 0)
#                xj = X[j] # position of current neighbour
#
#                # we need to find the cell index again, because this is
#                # not really the cell index, but it could be outside
#                # the domain -> i.e. this only makes a difference for pbc
#                cj = position_to_cell_index(inv_cell, xj, ns_vec)
#                cj = bin_trunc.(cj, pbc, ns_vec)
#
#                # drj is position relative to lower left corner of the bin
#                dxj = xj - bins * (cj - 1)
#                # Compute distance between atoms
#                dx = dxj - dxi + off
#                norm_dx_sq = dx ⋅ dx
#
#                # append to the list
#                if norm_dx_sq < cutoff_sq
#                   push!(first, i)
#                   push!(secnd, j)
#                   push!(distvec, dx)
#                   push!(absdist, sqrt(norm_dx_sq))
#                   push!(shift, (ci0 - cj + xyz) .÷ ns_vec)
#                end
#             end  # if i != j || any(xyz .!= 0)
#
#             # go to the next atom in the current cell
#             j = next[j];
#          end # while j > 0 (loop over atoms in current cell)
#       end # loop over neighbouring bins
#    end # for i = 1:nat
#
#    # Build return tuple
#    return first, secnd, absdist, distvec, shift
# end
