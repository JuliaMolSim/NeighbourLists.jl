
using NeighbourLists: SMat, SVec
using Test
using LinearAlgebra


include("nn_list.jl")

Ns = [3,3,4,4,5,5]
@info("Testing PairList Correctness: ")
for N in Ns
   C = SMat( diagm(0 => (2.0 .+ 0.2 * rand(3))) * N )
   X = [ C' * rand(SVec)   for i = 1:ceil(Int, abs(det(C))) ÷ 4 + 2 ]
   pbc = SVec(rand(Bool, 3))
   cutoff = 2.0

   # compute a cell list
   nlist = PairList(X, cutoff, C, pbc; sorted = true)

   # compute a NearestNeighbors list
   i, j, r, R = nn_list(X, cutoff, C, pbc)
   R = X[j] - X[i]
   first = NeighbourLists.get_first(i, length(X))
   NeighbourLists.sort_neigs!(j, r, R, zeros(SVec{Int32}, length(r)), first)

   println(@test (nlist.i == i) && (nlist.j == j))

   # check that they are sorted
   pass_sorted = true
   for (_1, j, _2, _3) in NeighbourLists.sites(nlist)
      if !issorted(j)
         pass_sorted = false
         break
      end
   end
   print("sorted: "); println(@test pass_sorted)

   # check that the neighbourhoods produced by `neigs` are correct
   pass_neigs = true
   for (i, j, r, R) in NeighbourLists.sites(nlist)
      j_, r_, R_ = neigs(nlist, i)
      if j != j_ || r != r_ || R != R_
         pass_neigs = false
         break
      end
   end
   print("neigs: "); println(@test pass_neigs)
end
println()
