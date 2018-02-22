using NeighbourList
using Base.Test
using JuLIP

function test_nlist(at, cutoff)
   cl = CellList(positions(at), cutoff, cell(at), pbc(at))
   nlist = neighbourlist(at, cutoff)
   @test (cl.i == nlist.i) && (cl.j == nlist.j) && (cl.S == nlist.S)
end

@testset "NeighbourList" begin

println("TEST 1: si, cubic cell, cluster")
at = bulk("Si", cubic=true) * 3
set_pbc!(at, false)
test_nlist(at, 1.1 * rnn("Si"))
test_nlist(at, 1.2 * rnn("Si"))

println("TEST 2: si, non-cubic cell, cluster")
at = bulk("Si") * 5
set_pbc!(at, false)
test_nlist(at, 2.1 * rnn("Si"))

println("TEST 3: si, non-cubic cell, pbc")
at = bulk("Si") * 5
set_pbc!(at, true)
test_nlist(at, 2.1 * rnn("Si"))

println("TEST 4: si, non-cubic cell, mixed bc")
at = bulk("Si") * 5
set_pbc!(at, (true, false, true))
test_nlist(at, 2.1 * rnn("Si"))

end

println("Performance Test:")
include("profile.jl")
