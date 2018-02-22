
import Base.mapreduce

export mapreduce!, mapreduce_d


function mapreduce!{T}(Es::Vector{T}, f, it::PairIterator{T})
   nlist = it.nlist
   @simd for n = 1:npairs(nlist)
      @inbounds Es[nlist.i[n]] += f(nlist.r[n], nlist.R[n])
   end
   return out
end

function mapreduce!{T}(out::Vector{SVec{T}}, df, it::PairIterator{T})
   nlist = it.nlist
   @simd for n = 1:npairs(nlist)
      df_ = df(nlist.r[n], nlist.R[n])
      out[nlist.j[n]] += df_
      out[nlist.i[n]] -= df_
   end
   return out
end

mapreduce{T}(f, it::PairIterator{T}) =
      mapreduce!(zeros(T, nsites(it.nlist)), f, it)

mapreduce_d{T}(df, it::PairIterator{T}) =
      mapreduce!(zeros(SVec{T}, nsites(it.nlist)), df, it)
