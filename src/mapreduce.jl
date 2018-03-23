
using Base.Threads

import Base.map!

export mapreduce!, mapreduce_sym!, mapreduce_antisym!, map_d!,
   mapreduce_d!, mapreduce_sym_d!

function mt_split(niter::TI, maxthreads=1_000_000_000) where TI
   nt = minimum([maxthreads, nthreads(), niter])
   nn = ceil.(TI, linspace(1, niter+1, nt+1))
   return nt, nn
end


"""
mapreduce!{S, T}(out::AbstractVector{S}, f, it::PairIterator{T})

basically a bin_sum, iterate over all pairs, for each pair evaluate
f(r, R) and add it to out[i] (the first)
"""
function mapreduce!{S, T}(out::AbstractVector{S}, f, it::PairIterator{T})
   nlist = it.nlist
   for n = 1:npairs(nlist)
      out[nlist.i[n]] += f(nlist.r[n], nlist.R[n])
   end
   return out
end


"""
`mapreduce_sym!{S, T}(out::AbstractVector{S}, f, it::PairIterator{T})`

symmetric variant of `mapreduce!{S, T}(out::AbstractVector{S}, ...)`, summing only
over bonds (i,j) with i < j and adding f(R_ij) to both sites i, j.
"""
function mapreduce_sym!{S, T}(out::AbstractVector{S}, f, it::PairIterator{T})
   nlist = it.nlist
   for n = 1:npairs(nlist)
      if  nlist.i[n] < nlist.j[n]    # NB this probably prevents simd
         f_ = f(nlist.r[n], nlist.R[n])
         out[nlist.i[n]] += f_
         out[nlist.j[n]] += f_
      end
   end
   return out
end

"""
`tmapreduce_sym!`

like `mapreduce_sym!` but multi-threaded
"""
function tmapreduce_sym!{S, T}(out::AbstractVector{S}, f, it::PairIterator{T})
   nlist = it.nlist
   nt, nn = mt_split(npairs(nlist))
   @threads for it = 1:nt
      for n = nn[it]:(nn[it+1]-1)
         if  nlist.i[n] < nlist.j[n]    # NB this probably prevents simd
            f_ = f(nlist.r[n], nlist.R[n])
            out[nlist.i[n]] += f_
            out[nlist.j[n]] += f_
         end
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
function mapreduce_antisym!{S, T}(out::AbstractVector{S}, f, it::PairIterator{T})
   nlist = it.nlist
   for n = 1:npairs(nlist)
      if nlist.i[n] < nlist.j[n]
         f_ = f(nlist.r[n], nlist.R[n])
         out[nlist.j[n]] += f_
         out[nlist.i[n]] -= f_
      end
   end
   return out
end

mapreduce_d!{S, T}(out::AbstractVector{S}, f, it::PairIterator{T}) =
   mapreduce_antisym!(out, f, it)


# ============ assembly over sites

function map!{S,T}(f, out::AbstractVector{S}, it::SiteIterator{T})
   nlist = it.nlist
   nt, nn = mt_split(nsites(nlist))
   @threads for it = 1:nt
      for i = nn[it]:(nn[it+1]-1)
         j, r, R = site(nlist, i)
         out[i] = f(r, R)
      end
   end
   return out
end

function map_d!{S,T}(f, out::AbstractVector{S}, it::SiteIterator{T})
   nlist = it.nlist
   nt, nn = mt_split(nsites(nlist))
   @threads for it = 1:nt
      for i = nn[it]:(nn[it+1]-1)
         j, r, R = site(nlist, i)
         f_ = f(r, R)
         out[j] += f_
         out[i] -= sum(f_)
      end
   end
   return out
end



# ============ assembly over n-body terms

"""
`@symm`: symmetrises a loop over a cartesian range. For example
```Julia
for i1 = a0:a1-2, i2 = i1+1:a1-1, i3 = i2+1:a1
   dosomething(i1, i2, i3)
end
```
may be written as
```Julia
@symm 3 for i = a0:a1
   dosomething(i[1], i[2], i[3])
end
```
here, `i` is a `CartesianIndex`.
"""
macro symm(NN, ex)
   N = eval(NN)
   @assert ex.head == :for
   @assert length(ex.args) == 2
   ex_for = ex.args[1]
   ex_body = ex.args[2]
   # iteration symbol
   i = ex_for.args[1]
   # lower and upper bound
   a0 = ex_for.args[2].args[1]
   a1 = ex_for.args[2].args[2]
   # create the for-loop without body
   loopstr = "for $(i)1 = ($a0):(($a1)-$(N-1))"
   for n = 2:N
      loopstr *= ", $i$n = $i$(n-1)+1:(($a1)-$(N-n))"
   end
   loopstr *= "\n $i = CartesianIndex($(i)1"
   for n = 2:N
      loopstr *= ", $i$n"
   end
   loopstr *= ") \n end"
   loopex = parse(loopstr)
   append!(loopex.args[2].args, ex_body.args)
   # return the expression
   esc(quote
      $loopex
   end)
end


struct NBodyIterator{N, T, TI}
   nlist::PairList{T,TI}
   order::Val{N}
end


"""
`nbodies(N, nlist::PairList)`

creates an N-body iterator, e.g., `nbodies(3, nlist)` will create an iterator
over 3-body terms. Use `mapreduce_sym!` and `map_reduce_sym_d!` to carry out the
iterations.
"""
nbodies(N, nlist::PairList) = NBodyIterator(nlist, Val(N))


"""
`function _find_next_(j, n, first)`

* `j` : array of neighbour indices
* `n` : current site index
* `first` : array of first indices

return the first index `first[n] <= m < first[n+1]` such that `j[m] > n`;
and returns 0 if no such index exists
"""
function _find_next_{TI}(j::Vector{TI}, n::TI, first::Vector{TI})
   # DEBUG CODE
   # @assert issorted(j[first[n]:first[n+1]-1])
   for m = first[n]:first[n+1]-1
      if j[m] > n
         return m
      end
   end
   return zero(TI)
end


function mapreduce_sym!(f, out::AbstractVector, it::NBodyIterator{N}) where N
   nlist = it.nlist
   nt, nn = mt_split(nsites(nlist))
   @threads for it = 1:nt
      for i = nn[it]:(nn[it+1]-1)
         # get the index of a neighbour > n
         a0 = _find_next_(nlist.j, n, nlist.first)
         a0 == 0 && continue  # (if no such index exists)
         # get the index up to which to loop
         a1 = nlist.first[n+1]-1
         j, r, R = site(nlist, i)
         @symm N for J = a0:a1
            # compute the N(N+1)/2 vector of distances
            s = simplex_lengths(r, R)
            out[j[J]] .+= f(s) / length(J)
         end
      end
   end
   return out
end


# function mapreduce_sym_d!{S,T}(f, out::AbstractVector{S}, it::NBodyIterator{3,T})
#    i, j, r, R, first = it.nlist.i, it.nlist.j, it.nlist.r, it.nlist.R, it.nlist.first
#    for n = 1:nsites(it.nlist)
#       # get the index of a neighbour > n
#       a0 = _find_next_(j, n, first)
#       if a0 == 0; continue; end  # (if no such index exists)
#       # get the index up to which to loop
#       a1 = first[n+1]-1
#       # loop over unique ordered tuples
#       for a = a0:a1-1, b = a0+1:a1
#          rab = norm(R[b]-R[a])
#          s = SVector{3, T}(r[a], rab, r[b])
#          f_ = f(s) / 3.0
#          f1a = (f_[1]/r[a]) * R[a]
#          f2ab = (f_[2]/rab) * (R[a]-R[b])
#          f3b = (f_[3]/r[b]) * R[b]
#          out[n] += - f1a  - f3b
#          out[j[a]] += f1a + f2ab
#          out[j[b]] += - f2ab + f3b
#       end
#    end
#    return out
# end
