using NeighbourList
using Base.Test
using JuLIP

# performance = false
#
# function test_nlist(at, cutoff)
#    cl = CellList(positions(at), cutoff, cell(at), pbc(at))
#    nlist = neighbourlist(at, cutoff)
#    @test (cl.i == nlist.i) && (cl.j == nlist.j) && (cl.S == nlist.S)
# end
#
#
# @testset "NeighbourList" begin
#    include("test_celllist.jl")
# end
#
#
# if performance
#    println("`NeighbourLists` Performance Tests:")
#    include("profile.jl")
# end


include("test_forces.jl")
