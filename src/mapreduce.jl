
using Base.Threads

import Base.map!

export mapreduce!, mapreduce_sym!, mapreduce_antisym!, map_d!,
   mapreduce_d!, mapreduce_sym_d!

function mt_split(niter::TI, maxthreads=1_000_000_000) where TI
   nt = minimum([maxthreads, nthreads(), niter])
   nn = ceil.(TI, linspace(1, niter+1, nt+1))
   return nt, nn
end

function mt_split_interlaced(niter::TI, maxthreads=1_000_000_000) where TI
   nt = minimum([maxthreads, nthreads(), niter])
   nn = [ j:nt:niter for j = 1:nt ]
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
macro symm(N, ex)
   if N isa Symbol
      N = eval(N)
   end
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
   loopstr *= "\n $i = SVector{$N, Int}($(i)1"
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

"""
`simplex_lengths`: compute the sidelengths of a simplex
and return the corresponding pairs of X indices
"""
function simplex_lengths!(s, a, b, i, J::SVector{N, TI}, nlist
                           ) where {N, TI <: Integer}
   n = 0
   for l = 1:N
      n += 1
      a[n] = i
      b[n] = nlist.j[J[l]]
      s[n] = nlist.r[J[l]]
   end
   for i1 = 1:N-1, j1 = i1+1:N
      n += 1
      a[n] = nlist.j[J[i1]]
      b[n] = nlist.j[J[j1]]
      s[n] = norm(nlist.R[J[i1]] - nlist.R[J[j1]])
   end
   return SVector(s), SVector(a), SVector(b)
end


@generated function mr_sym_inner!(f, out::AbstractVector,
         it::NBodyIterator{N, T, TI}, rg) where {N, T, TI}
   N2 = (N*(N-1))÷2
   quote
      nlist = it.nlist
      # allocate some temporary arrays
      a_ = zero(MVector{$N2, TI})
      b_ = zero(MVector{$N2, TI})
      s_ = zero(MVector{$N2, T})
      # loop over the range allocated to this thread
      for i in rg
         # get the index of a neighbour > n
         a0 = _find_next_(nlist.j, i, nlist.first)
         a0 == 0 && continue  # (if no such index exists)
         # get the index up to which to loop
         a1 = nlist.first[i+1]-1
         @symm $(N-1) for J = a0:a1
            # compute the N(N+1)/2 vector of distances
            s, _, _ = simplex_lengths!(s_, a_, b_, i, J, nlist)
         #                        ~~~~~~~~~~~~~~~~~~~ generic up to here
            f_ = f(s) / $N
            out[i] += f_
            for l = 1:length(J)
               out[nlist.j[J[l]]] .+= f_
            end
         end
      end
   end
end

function mapreduce_sym!(
         f, out::AbstractVector, it::NBodyIterator{N, T, TI}) where {N, T, TI}
   nt, nn = mt_split(nsites(it.nlist))
   for i = 1:nt
      mr_sym_inner!(f, out, it, nn[i]:(nn[i+1]-1))
   end
   return out
end


@generated function mr_sym_d_inner!(df, out::AbstractVector,
         it::NBodyIterator{N, T, TI}, rg) where {N, T, TI}
   N2 = (N*(N-1))÷2
   quote
      nlist = it.nlist
      # allocate some temporary arrays
      a_ = zero(MVector{$N2, TI})
      b_ = zero(MVector{$N2, TI})
      s_ = zero(MVector{$N2, T})
      # loop over the range allocated to this thread
      for i in rg
         # get the index of a neighbour > n
         a0 = _find_next_(nlist.j, i, nlist.first)
         a0 == 0 && continue  # (if no such index exists)
         # get the index up to which to loop
         a1 = nlist.first[i+1]-1
         @symm $(N-1) for J = a0:a1
            # compute the N(N+1)/2 vector of distances
            s, a, b = simplex_lengths!(s_, a_, b_, i, J, nlist)
         #                        ~~~~~~~~~~~~~~~~~~~ generic up to here
            df_ = df(s)
            for l = 1:length(s)
               Rab = nlist.X[a[l]] - nlist.X[b[l]]
               Sab = Rab / norm(Rab)
               out[a[l]] += df_[l] * Sab
               out[b[l]] -= df_[l] * Sab
            end
         end
      end
   end
end

function mapreduce_sym_d!(
         df, out::AbstractVector, it::NBodyIterator{N, T, TI}) where {N, T, TI}
   nt, nn = mt_split(nsites(it.nlist))
   for i = 1:nt
      mr_sym_d_inner!(df, out, it, nn[i]:(nn[i+1]-1))
   end
   return out
end
