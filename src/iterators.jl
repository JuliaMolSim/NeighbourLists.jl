
import Base: start, done, next, length

export pairs, sites, site, nbodies

inc{T <: Integer}(i::T) = i + T(1)

# -------------- iterator over pairs ---------------

pairs(nlist::PairList) = PairIterator(nlist)


struct PairIterator{T,TI}
   nlist::PairList{T,TI}
end

start{T,TI}(it::PairIterator{T,TI}) = TI(1)
done(it::PairIterator, i::Integer) = (i > npairs(it.nlist))
next(it::PairIterator, i) =
   (it.nlist.i[i], it.nlist.j[i], it.nlist.r[i], it.nlist.R[i]), inc(i)
length(it::PairIterator) = length(it.nlist)

# -------------- iterator over sites ---------------

function site(nlist::PairList, i0)
   n1, n2 = nlist.first[i0], nlist.first[i0+1]-1
   return (@view nlist.j[n1:n2]), (@view nlist.r[n1:n2]), (@view nlist.R[n1:n2])
end

sites(nlist::PairList) = SiteIterator(nlist)

struct SiteIterator{T,TI}
   nlist::PairList{T,TI}
end

start{T,TI}(it::SiteIterator{T,TI}) = one(TI)
done(it::SiteIterator, i::Integer) = (i > nsites(it.nlist))
next(it::SiteIterator, i::Integer) = (i, site(it.nlist, i)...), inc(i)
length(it::SiteIterator) = nsites(it.nlist)


# -------------- iterator over n-body terms ---------------


struct NBodyIterator{N, T, TI}
   nlist::PairList{T,TI}
   order::Val{N}
end


"""
`nbodies(N, nlist::PairList)`

creates an N-body iterator, e.g., `nbodies(3, nlist)` will create an iterator
over permutation-symmetric 3-body terms. Use `mapreduce_sym!`
and `map_reduce_sym_d!` to carry out the iterations.
"""
nbodies(N::Integer, nlist::PairList) = NBodyIterator(nlist, Val(N))
