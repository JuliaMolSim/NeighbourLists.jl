using NeighbourLists: SMat, SVec
Ns = [2,2,3,3,4,4,5,5]
for N in Ns
   C = SMat( (2 * eye(3) + 0.5 * sign.(rand(3,3)-0.5) + 0.1 * rand(3,3)-0.1) ) * N
   X = [ C * rand(SVec)   for i = 1:ceil(Int, abs(det(C))) + 2 ]
   @show det(C), length(X)
   pbc = (false, false, false)
   cutoff = 2.0
   nlist = CellList(X, cutoff, X, pbc)
   nnlist =
end

println("testing sort_neigs!")
_, at, cutoff = test_configs[2]
nlist = CellList(positions(at), cutoff, cell(at), pbc(at), sorted=true)
pass_sorted = true
for (_1, j, _2, _3) in NeighbourLists.sites(nlist)
   if !issorted(j)
      pass_sorted = false
      break
   end
end
@test pass_sorted


using NeighbourLists: SMat
sign.(rand(3,3) - 0.5) +


abs.(det.(cells))
