
import Base: start, done, next

export pairs, sites, site


# -------------- iterator over pairs ---------------

pairs(nlist::CellList) = PairIterator(nlist)

struct PairIterator{T,TI}
   nlist::CellList{T,TI}
end
start{T,TI}(it::PairIterator{T,TI}) = one(TI)
done(it::PairIterator, i::Integer) = (i == npairs(it.nlist) + 1)
next{T,TI}(it::PairIterator{T,TI}, i::TI) =
   it.nlist.i[i], it.nlist.j[i], it.nlist.r[i], it.nlist.R[i]

# -------------- iterator over sites ---------------

function site(nlist::CellList, i0)
   n1, n2 = nlist.first[i0], nlist.first[i0+1]-1
   return @views (nlist.j[n1:n2], nlist.r[n1:n2], nlist.R[n1:n2]) 
end

sites(nlist::CellList) = SiteIterator(nlist)

struct SiteIterator{T,TI}
   nlist::CellList{T,TI}
end


# -------------- iterator over bond-angles ---------------
