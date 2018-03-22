
using Base.Threads

import Base.map!

export mapreduce!, mapreduce_sym!, mapreduce_antisym!, map_d!,
   mapreduce_d!, mapreduce_sym_d!


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
   nperthread = npairs(nlist) ÷ nthreads()
   nn = [1:nperthread:npairs(nlist); npairs(nlist)+1]
   @threads for it = 1:nthreads()
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
   for i = 1:nsites(nlist)
      j, r, R = site(nlist, i)
      out[i] = f(r, R)
   end
   return out
end

function map_d!{S,T}(f, out::AbstractVector{S}, it::SiteIterator{T})
   nlist = it.nlist
   for i = 1:nsites(nlist)
      j, r, R = site(nlist, i)
      f_ = f(r, R)
      out[j] += f_
      out[i] -= sum(f_)
   end
   return out
end



# ============ assembly over n-body terms

struct NBodyIterator{N, T, TI}
   nlist::CellList{T,TI}
   order::Val{N}
end


"""
`function _find_next_(j, n, first)`

* `j` : array of neighbour indices
* `n` : current site index
* `first` : array of first indices

return the first index `first[n] <= m < first[n]` such that `j[m] > n`;
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

function mapreduce_sym!{S,T}(f, out::AbstractVector{S}, it::NBodyIterator{3,T})
   nlist = it.nlist
   for n = 1:nsites(nlist)
      # get the index of a neighbour > n
      a0 = _find_next_(nlist.j, n, nlist.first)
      if a0 == 0; continue; end  # (if no such index exists)
      # get the index up to which to loop
      a1 = nlist.first[n+1]-1
      # loop over unique ordered tuples
      for a = a0:a1-1, b = a0+1:a1
         s = SVector{3, T}(r[a], norm(R[b]-R[a]), r[b])
         f_ = f(s) / 3.0
         out[n] += f_
         out[nlist.j[a]] += f_
         out[nlist.j[b]] += f_
      end
   end
   return out
end


function mapreduce_sym_d!{S,T}(f, out::AbstractVector{S}, it::NBodyIterator{3,T})
   i, j, r, R, first = it.nlist.i, it.nlist.j, it.nlist.r, it.nlist.R, it.nlist.first
   for n = 1:nsites(it.nlist)
      # get the index of a neighbour > n
      a0 = _find_next_(j, n, first)
      if a0 == 0; continue; end  # (if no such index exists)
      # get the index up to which to loop
      a1 = first[n+1]-1
      # loop over unique ordered tuples
      for a = a0:a1-1, b = a0+1:a1
         rab = norm(R[b]-R[a])
         s = SVector{3, T}(r[a], rab, r[b])
         f_ = f(s) / 3.0
         f1a = (f_[1]/r[a]) * R[a]
         f2ab = (f_[2]/rab) * (R[a]-R[b])
         f3b = (f_[3]/r[b]) * R[b]
         out[n] += - f1a  - f3b
         out[j[a]] += f1a + f2ab
         out[j[b]] += - f2ab + f3b
      end
   end
   return out
end
