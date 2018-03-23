

using NeighbourLists
using JuLIP


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


function s(b0, b1)
   @symm 3 for j = b0:b1
      println(j)
   end
end

s(3, 8)


function t(b0, b1, ::Val{N}) where N
   @symm N for j = b0:b1
      println(j)
   end
end

t(3, 8, Val(3))
