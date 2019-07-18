
import Base: iterate, length, pairs

export pairs, sites, site, nbodies

abstract type AbstractIterator end

inc(i::T) where {T <: Integer} = i + T(1)

# -------------- iterator over pairs ---------------

pairs(nlist::PairList) = PairIterator(nlist)


struct PairIterator{T,TI} <: AbstractIterator
   nlist::PairList{T,TI}
end

_item(it::PairIterator, i::Integer) =
      (it.nlist.i[i], it.nlist.j[i], _getR(it.nlist, i))
iterate(it::PairIterator{T,TI}) where {T,TI} =
   npairs(it.nlist) > 0 ? (_item(it, 1), TI(1)) : nothing
iterate(it::PairIterator, i::Integer) =
   i >= npairs(it.nlist) ? nothing : (_item(it, inc(i)), inc(i))
length(it::PairIterator) = npairs(it.nlist)

# -------------- iterator over sites ---------------


sites(nlist::PairList) = SiteIterator(nlist)

struct SiteIterator{T,TI}  <: AbstractIterator
   nlist::PairList{T,TI}
end

_item(it::SiteIterator, i::Integer) = (i, neigs(it.nlist, i)...)
iterate(it::SiteIterator{T,TI}) where {T,TI} = _item(it, 1), one(TI)
iterate(it::SiteIterator, i::Integer) =
   i >= length(it) ? nothing : (_item(it, i+1), inc(i))
length(it::SiteIterator) = nsites(it.nlist)
