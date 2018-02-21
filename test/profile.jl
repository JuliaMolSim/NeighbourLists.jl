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

# si, non-cubic cell, mixed bc
at = bulk("Si", cubic=true) * 10
set_pbc!(at, (true, false, true))
C = JMat(cell(at))
X = positions(at)
perbc = JVec(pbc(at))
cutoff = 2.1 * rnn("Si")

println("Julia Nlist")
@btime NeighbourList.neighbour_list(C, perbc, X, cutoff)
println("Matscipy Nlist")
@btime matscipy_nlist(at, cutoff)
