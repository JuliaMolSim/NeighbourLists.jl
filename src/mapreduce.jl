
using Base.Threads

export maptosites!, maptosites_d!


function mt_split(niter::TI, maxthreads=MAX_THREADS[1]) where TI
   nt = minimum([maxthreads, nthreads(), niter])
   # nn = ceil.(TI, linspace(1, niter+1, nt+1))
   nn = ceil.(TI, range(1, stop=niter+1, length=nt+1))
   rgs = [nn[i]:(nn[i+1]-1) for i = 1:nt]
   return nt, rgs
end

function mt_split_interlaced(niter::TI, maxthreads=MAX_THREADS[1]) where TI
   nt = minimum([maxthreads, nthreads(), niter])
   rgs = [ j:nt:niter for j = 1:nt ]
   return nt, rgs
end


function _mt_map_!(f::FT, out, it, inner_loop) where FT
   nt, rg = mt_split(length(it))
   if nt == 1
      inner_loop(f, out, it, 1:length(it))
   else
      OUT = [[out]; [zeros(out) for i = 2:nt]]
      @threads for i = 1:nt
         inner_loop(f, OUT[i], it, rg[i])
      end
      for it = 2:nt
         out .+= OUT[it]
      end
   end
   return out
end

maptosites!(f, out::AbstractVector, it::AbstractIterator) =
   _mt_map_!(f, out, it, maptosites_inner!)

maptosites_d!(f, out::AbstractVector, it::AbstractIterator) =
   _mt_map_!(f, out, it, maptosites_d_inner!)

# ============ assembly over pairs


"""
`mapreduce_sym!{S, T}(out::AbstractVector{S}, f, it::PairIterator{T})`

symmetric variant of `mapreduce!{S, T}(out::AbstractVector{S}, ...)`, summing only
over bonds (i,j) with i < j and adding f(R_ij) to both sites i, j.
"""
function maptosites_inner!(f::FT, out, it::PairIterator, rg) where FT
   nlist = it.nlist
   for n in rg
      if  nlist.i[n] < nlist.j[n]
         f_ = f(nlist.r[n], nlist.R[n]) / 2
         out[nlist.i[n]] += f_
         out[nlist.j[n]] += f_
      end
   end
   return out
end


"""
`mapreduce_antisym!{T}(out::AbstractVector{SVec{T}}, df, it::PairIterator{T})`

anti-symmetric variant of `mapreduce!{S, T}(out::AbstractVector{S}, ...)`, summing only
over bonds (i,j) with i < j and adding f(R_ij) to site j and
-f(R_ij) to site i.
"""
function maptosites_d_inner!(f::FT, out, it::PairIterator, rg) where FT
   nlist = it.nlist
   for n in rg
      if nlist.i[n] < nlist.j[n]
         f_ = f(nlist.r[n], nlist.R[n])
         out[nlist.j[n]] += f_
         out[nlist.i[n]] -= f_
      end
   end
   return out
end


# ============ assembly over sites

function maptosites!(f::FT, out::AbstractVector, it::SiteIterator) where FT
   @threads for i = 1:nsites(it.nlist)
      _, r, R = neigs(it.nlist, i)
      out[i] = f(r, R)
   end
   return out
end

function maptosites_d!(df::FT, out::AbstractVector, it::SiteIterator) where FT
   nt = nthreads()
   OUT = [out; [zeros(out) for n = 2:nt]]
   @threads for i = 1:nsites(it.nlist)
      j, r, R = neigs(it.nlist, i)
      df_ = df(r, R)
      OUT[threadid()][j] += df_
      OUT[threadid()][i] -= sum(df_)
   end
   for it = 2:n
      out .+= OUT[it]
   end
   return out
end



# # ============ assembly over n-body terms
#
# """
# `@symm`: symmetrises a loop over a cartesian range. For example
# ```Julia
# for i1 = a0:a1-2, i2 = i1+1:a1-1, i3 = i2+1:a1
#    dosomething(i1, i2, i3)
# end
# ```
# may be written as
# ```Julia
# @symm 3 for i = a0:a1
#    dosomething(i[1], i[2], i[3])
# end
# ```
# here, `i` is a `CartesianIndex`.
# """
# macro symm(N, ex)
#    if N isa Symbol
#       N = eval(N)
#    end
#    @assert ex.head == :for
#    @assert length(ex.args) == 2
#    ex_for = ex.args[1]
#    ex_body = ex.args[2]
#    # iteration symbol
#    i = ex_for.args[1]
#    # lower and upper bound
#    a0 = ex_for.args[2].args[1]
#    a1 = ex_for.args[2].args[2]
#    # create the for-loop without body, e.g., for a 3-body assembly it generates
#    #   for i1 = a0:a1-2, i2 = i1+1:a1-1, i3 = i2+1:a1
#    #      do something with (i1, i2, i3)
#    #   end
#    loopstr = "for $(i)1 = ($a0):(($a1)-$(N-1))"
#    for n = 2:N
#       loopstr *= ", $i$n = $i$(n-1)+1:(($a1)-$(N-n))"
#    end
#    loopstr *= "\n $i = SVector{$N, Int}($(i)1"
#    for n = 2:N
#       loopstr *= ", $i$n"
#    end
#    loopstr *= ") \n end"
#    loopex = Meta.parse(loopstr)
#    append!(loopex.args[2].args, ex_body.args)
#    # return the expression
#    esc(quote
#       $loopex
#    end)
# end




"""
`function _find_next_(j, n, first)`

* `j` : array of neighbour indices
* `n` : current site index
* `first` : array of first indices

return the first index `first[n] <= m < first[n+1]` such that `j[m] > n`;
and returns 0 if no such index exists
"""
function _find_next_(j::Vector{TI}, n::TI, first::Vector{TI}) where TI
   # DEBUG CODE
   # @assert issorted(j[first[n]:first[n+1]-1])
   for m = first[n]:first[n+1]-1
      if j[m] > n
         return m
      end
   end
   return zero(TI)
end

# """
# `simplex_lengths`: compute the sidelengths of a simplex
# and return the corresponding pairs of X indices
# """
# function simplex_lengths!(s, S, a, b, i, J::SVector{N, TI}, nlist
#                            ) where {N, TI <: Integer}
#    n = 0
#    for l = 1:N
#       n += 1
#       a[n] = i
#       b[n] = nlist.j[J[l]]
#       s[n] = nlist.r[J[l]]
#       S[n] = - nlist.R[J[l]] / s[n]  # Rab
#    end
#    for i1 = 1:N-1, j1 = i1+1:N
#       n += 1
#       a[n] = nlist.j[J[i1]]
#       b[n] = nlist.j[J[j1]]
#       Rab = nlist.R[J[i1]] - nlist.R[J[j1]]
#       s[n] = norm(Rab)
#       S[n] = Rab / s[n]
#    end
#    return SVector(s), SVector(S), SVector(a), SVector(b)
# end

# _m2s_mul_(x::T, S::SVector{N,T}) where {N, T} = x * S

# @inline function _inc_stress_!(out::MMatrix, s, df::Number, S)
#    out[:, :] -= (s * df) * (S * S')
# end


# @generated function _m2s_generic_!(f::FT, out::AbstractArray,
#                         it::NBodyIterator{N, T, TI}, rg,
#                         forgrad::Val{FG}) where {FT, N, T, TI, FG}
#    N2 = (N*(N-1))÷2
#    if FG == :F   # energy
#       mapcode = quote
#          f_ = f(s) / $N
#          out[i] += f_
#          for l = 1:length(J)
#             out[nlist.j[J[l]]] += f_
#          end
#       end
#    elseif FG == :G  # gradient
#       mapcode = quote
#          df_ = f(s)
#          for l = 1:length(s)
#             _t = _m2s_mul_(df_[l], S[l])   # df_[l] * S[l]
#             out[a[l]] += _t
#             out[b[l]] -= _t
#          end
#       end
#    elseif FG == :V   # virial (f = ∇Vn(r12, ...))
#       mapcode = quote
#          df_ = f(s)
#          for l = 1:length(s)
#             # out[:, :] -= (s[l]*df_[l]) * (S[l] * S[l]')
#             _inc_stress_!(out, s[l], df_[l], S[l])
#          end
#       end
#    else
#       error("unknown `FG`")
#    end
#
#    quote
#       nlist = it.nlist
#       # allocate some temporary arrays
#       a_ = zero(MVector{$N2, TI})
#       b_ = zero(MVector{$N2, TI})
#       s_ = zero(MVector{$N2, T})
#       S_ = zero(MVector{$N2, SVec{T}})
#       # loop over the range allocated to this thread
#       for i in rg
#          # get the index of a neighbour > n
#          a0 = _find_next_(nlist.j, i, nlist.first)
#          a0 == 0 && continue  # (if no such index exists)
#          # get the index up to which to loop
#          a1 = nlist.first[i+1]-1
#          @symm $(N-1) for J = a0:a1
#             # compute the N(N+1)/2 vector of distances
#             s, S, a, b = simplex_lengths!(s_, S_, a_, b_, i, J, nlist)
#             if maximum(s) < nlist.cutoff
#                $mapcode
#             end
#          end
#       end
#    end
# end


# maptosites_inner!(f, out, it::NBodyIterator, rg) =
#       _m2s_generic_!(f, out, it, rg, Val(:F))
#
# maptosites_d_inner!(f, out, it::NBodyIterator, rg) =
#       _m2s_generic_!(f, out, it, rg, Val(:G))
#
# virial!(f, out, it::NBodyIterator) =
#       _m2s_generic_!(f, out, it, 1:nsites(it.nlist), Val(:V))
