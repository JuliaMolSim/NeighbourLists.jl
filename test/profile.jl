using NeighbourList
using JuLIP
using BenchmarkTools
using PyCall
using ProfileView

# @pyimport matscipy.neighbours as matscipy_neighbours
matscipy_neighbours = pyimport("matscipy.neighbours")

function matscipy_nlist(at, cutoff)
   pycall(matscipy_neighbours["neighbour_list"],
          NTuple{5, PyArray}, "ijdDS", at.po, cutoff)
end

# TEST 4: si, non-cubic cell, mixed bc
at = bulk("Si", cubic=true) * 10
set_pbc!(at, (true, false, true))
C = JMat(cell(at))
X = positions(at)
perbc = JVec(pbc(at))
cutoff = 2.1 * rnn("Si")

# NeighbourList.neighbour_list(C, perbc, X, cutoff)

println("Julia Nlist-v1")
@btime NeighbourList.neighbour_list(C, perbc, X, cutoff)
# println("Julia Nlist-v2")
# @btime NeighbourList.neighbour_list2(C, perbc, X, cutoff, Int)
println("Matscipy Nlist")
@btime matscipy_nlist(at, cutoff)



# =================== TEMP


# NeighbourList.neighbour_list2(C, perbc, X, cutoff);
#
# Base.Profile.clear()
# @profile NeighbourList._neighbour_list2_(C, perbc, X, cutoff, 0);
#
# ProfileView.view()


# @btime neighbourlist(rattle!(at, 0.0001), cutoff)

# Base.Profile.clear()
# @profile NeighbourList.neighbour_list(C, perbc, X, cutoff);
# Base.Profile.print()

# @code_warntype NeighbourList.neighbour_list(C, perbc, X, cutoff)
