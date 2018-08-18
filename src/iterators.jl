
import Base: iterate, length

export pairs, sites, site, nbodies

abstract type AbstractIterator end

inc(i::T) where {T <: Integer} = i + T(1)

# -------------- iterator over pairs ---------------

pairs(nlist::PairList) = PairIterator(nlist)


struct PairIterator{T,TI} <: AbstractIterator
   nlist::PairList{T,TI}
end

_item(it::PairIterator, i::Integer) =
               (it.nlist.i[i], it.nlist.j[i], it.nlist.r[i], it.nlist.R[i])
iterate(it::PairIterator{T,TI}) where {T,TI} = _item(it, 1), TI(1)
iterate(it::PairIterator, i::Integer) =
   i >= npairs(it.nlist) ? nothing : (_item(it, inc(i)), inc(i))
length(it::PairIterator) = npairs(it.nlist)

# -------------- iterator over sites ---------------

function site(nlist::PairList, i0)
   n1, n2 = nlist.first[i0], nlist.first[i0+1]-1
   return (@view nlist.j[n1:n2]), (@view nlist.r[n1:n2]), (@view nlist.R[n1:n2])
end

sites(nlist::PairList) = SiteIterator(nlist)

struct SiteIterator{T,TI}  <: AbstractIterator
   nlist::PairList{T,TI}
end

_item(it::SiteIterator, i::Integer) = (i, site(it.nlist, i)...)
iterate(it::SiteIterator{T,TI}) where {T,TI} = _item(it, 1), one(TI)
iterate(it::SiteIterator, i::Integer) =
   i >= length(it) ? nothing : (_item(it, i+1), inc(i))
length(it::SiteIterator) = nsites(it.nlist)


# -------------- iterator over n-body terms ---------------


struct NBodyIterator{N, T, TI}  <: AbstractIterator
   nlist::PairList{T,TI}
   order::Val{N}
end

length(it::NBodyIterator) = nsites(it.nlist)

"""
`nbodies(N, nlist::PairList)`

creates an N-body iterator, e.g., `nbodies(3, nlist)` will create an iterator
over permutation-symmetric 3-body terms. Use `mapreduce_sym!`
and `map_reduce_sym_d!` to carry out the iterations.
"""
nbodies(N::Integer, nlist::PairList) = NBodyIterator(nlist, Val(N))
