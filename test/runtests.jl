using NeighbourList
using Base.Test
using JuLIP

function test_nlist(at, cutoff)
   C = JMat(cell(at))
   X = positions(at)
   perbc = JVec(pbc(at))
   i, j, r, R, S = NeighbourList.neighbour_list(C, perbc, X, cutoff)
   nlist = neighbourlist(at, cutoff)
   @test (r == nlist.r) && (R == nlist.R) && (S == nlist.S) &&
         (i == nlist.i) && (j == nlist.j)
end

@testset "NeighbourList" begin

# TEST 1: si, cubic cell, cluster
at = bulk("Si", cubic=true) * 3
set_pbc!(at, false)
test_nlist(at, 1.1 * rnn("Si"))
test_nlist(at, 1.2 * rnn("Si"))

# TEST 2: si, non-cubic cell, cluster
at = bulk("Si") * 5
set_pbc!(at, false)
test_nlist(at, 2.1 * rnn("Si"))

# TEST 3: si, non-cubic cell, pbc
at = bulk("Si") * 5
set_pbc!(at, true)
test_nlist(at, 2.1 * rnn("Si"))

# TEST 4: si, non-cubic cell, mixed bc
at = bulk("Si") * 5
set_pbc!(at, (true, false, true))
test_nlist(at, 2.1 * rnn("Si"))

end
