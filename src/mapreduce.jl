
import Base.map!

export mapreduce!, mapreduce_sym!, mapreduce_antisym!, map_cfd!


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


# ============ assembly over sites

function map!{S,T}(f, out::AbstractVector{S}, it::SiteIterator{T})
   nlist = it.nlist
   for i = 1:nsites(nlist)
      j, r, R = site(nlist, i)
      out[i] = f(r, R)
   end
   return out
end

function map_cfd!{S,T}(f, out::AbstractVector{S}, it::SiteIterator{T})
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

# function mapreduce_sym!{S,T}(f, out::AbstractVector{S}, it::NBodyIterator{3,T})
#    nlist = it.nlist
#    for (i, j, r, R) in sites(nlist)
#       for a = 1:length(j), b = 1:length(j)
#          if !(i < j[a] < j[b])
#             continue
#          end
#          s = SVector{3, T}(r[a], norm(R[b]-R[a]), r[b])
#          f_ = f(s) / 3.0
#          out[i] += f_
#          out[j[a]] += f_
#          out[j[b]] += f_
#       end
#    end
#    return out
# end
