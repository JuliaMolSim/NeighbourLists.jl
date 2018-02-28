
include("test_configs.jl")

@testset "CellList" begin

for (i, (descr, at, cutoff)) in enumerate(test_configs)
   println("TEST $i: $descr")
   test_nlist(at, cutoff)
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

end
