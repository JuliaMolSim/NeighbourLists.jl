

using NearestNeighbors, Distances, LinearAlgebra

import Distances: evaluate


struct Periodic{T} <: Metric
   L::SVec{T}
   pbc::SVec{Bool}
end

function evaluate(d::Periodic{T}, x::AbstractVector, y::AbstractVector) where {T}
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

# function nn_list{T, TI}(X::Vector{SVec{T}}, cutoff::T,
#                         cell::SMat{T}, pbc::SVec{Bool}, _::TI)
function nn_list(X::Vector{SVec{Float64}}, cutoff, cell, pbc::SVec{Bool})
   # WARNING: removed `const` here - will this slow down the code?
   T = Float64
   TI = Int
   @assert cell == SMat{T}(diagm(0 => diag(cell)))   # require cubic cell

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
            push!(r, evaluate(d, X[n], X[neigs[a]]))
         end
      end
   end

   return i, j, r, R
end
