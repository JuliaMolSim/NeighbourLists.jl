
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
      _, R = neigs(it.nlist, i)
      out[i] = f(R)
   end
   return out
end

function maptosites_d!(df::FT, out::AbstractVector, it::SiteIterator) where FT
   nt = nthreads()
   OUT = [out; [zeros(out) for n = 2:nt]]
   @threads for i = 1:nsites(it.nlist)
      j, R = neigs(it.nlist, i)
      df_ = df(R)
      OUT[threadid()][j] += df_
      OUT[threadid()][i] -= sum(df_)
   end
   for it = 2:n
      out .+= OUT[it]
   end
   return out
end
