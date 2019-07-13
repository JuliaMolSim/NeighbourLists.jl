using NeighbourLists
using Test

# ---- FLAGS -----

# whether to run performance tests
performance = true

# check whether on CI
isCI = haskey(ENV, "CI")
notCI = !isCI

# TODO: switch the JuLIP test to an ASE test
# check whether we have JuLIP
# hasjulip = true
# try
#    using JuLIP
# catch
#    hasjulip = false
# end

# ----------------- TESTS -------------------

println("# threads = $(Base.Threads.nthreads())")

@testset "NeighbourLists" begin
   @testset "Aux" begin include("test_aux.jl") end
   @testset "CellList" begin include("test_celllist.jl") end

   # pointless until we switch to comparing against ASE / matscipy
   # if hasjulip
   #    @testset "JuLIP" begin include("test_julip.jl") end
   # end
end
