
include("test_configs.jl")

@testset "CellList" begin

for (i, (descr, at, cutoff)) in enumerate(test_configs)
   println("TEST $i: $descr")
   test_nlist(at, cutoff)
end

end
