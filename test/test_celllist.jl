
using NeighbourLists: SMat, SVec
using Base.Test

include("nn_list.jl")

Ns = [3,3,4,4,5,5]
for N in Ns

   #  * eye(3) + 0.5 * sign.(rand(3,3)-0.5) + 0.1 * rand(3,3)-0.1) ) * N
   C = SMat( diagm(2.0 + 0.2 * rand(3)) * N )
   X = [ C' * rand(SVec)   for i = 1:ceil(Int, abs(det(C))) + 2 ]
   pbc = SVec(rand(Bool, 3))
   cutoff = 2.0

   # compute a cell list
   nlist = PairList(X, cutoff, C, pbc; sorted = true)

   # compute a NearestNeighbors list
   i, j, r, R = nn_list(X, cutoff, C, pbc)
   R = X[j] - X[i]
   first = NeighbourLists.get_first(i, length(X))
   NeighbourLists.sort_neigs!(j, r, R, first)

   # check the two lists are identical
   @test (i == nlist.i) && (j == nlist.j)

   # check that they are sorted
   pass_sorted = true
   for (_1, j, _2, _3) in NeighbourLists.sites(nlist)
      if !issorted(j)
         pass_sorted = false
         break
      end
   end
   @test pass_sorted
end
