

using NeighbourLists
using JuLIP
using Base.Threads

# at = set_pbc!(bulk(:Si, cubic=true) * 3, true)
# length(at)
# cutoff = 3.1 * rnn(:Si)
# nlist = CellList(positions(at), cutoff, cell(at), pbc(at))
#
# f(r, R) = (exp(-8*(r-1)) - 2 * exp( - 4 * (r - 1))) * 1 / (1+r)
#
# out = zeros(length(at))
# mapreduce_sym!(out, f, pairs(nlist))
#
# @time mapreduce_sym!(out, f, pairs(nlist));
# @time NeighbourLists.tmapreduce_sym!(out, f, pairs(nlist));



struct NBodyIterator{N, T, TI}
   nlist::PairList{T,TI}
   order::Val{N}
end

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


function t(it::Val{N})
   @symm N for i = 1:5
      println(i)
   end
end
