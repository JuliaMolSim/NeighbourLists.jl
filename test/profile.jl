using NeighbourLists
using JuLIP
using BenchmarkTools
using PyCall
using ASE
# using ProfileView

# @pyimport matscipy.neighbours as matscipy_neighbours
matscipy_neighbours = pyimport("matscipy.neighbours")

function matscipy_nlist(at, cutoff)
   pycall(matscipy_neighbours["neighbour_list"],
          NTuple{5, PyArray}, "ijdDS", at.po, cutoff)
end

println("# Threads = ", Base.Threads.nthreads())

for L in [4, 10, 30]
   print("L = $L")
   # si, non-cubic cell, mixed bc
   at = bulk(:Si, cubic=true) * L
   println(", N = $(length(at))")
   set_pbc!(at, (true, false, true))
   C = JMat(cell(at))
   X = positions(at)
   perbc = JVec(pbc(at))
   cutoff = 2.1 * rnn(:Si)

   println("Julia Nlist")
   @btime PairList($X, $cutoff, $C, $perbc, sorted = false)
   println("Matscipy Nlist")
   @btime matscipy_nlist($(ASEAtoms(at)), $cutoff)
   println("------------------------------------------")
end

using NeighbourLists
using JuLIP
using BenchmarkTools

L = 30
print("L = $L")
# si, non-cubic cell, mixed bc
at = bulk(:Si, cubic=true) * L
println(", N = $(length(at))")
set_pbc!(at, (true, false, true))
C = JMat(cell(at))
X = positions(at)
perbc = JVec(pbc(at))
cutoff = 2.1 * rnn(:Si)
nlist = PairList(X, cutoff, C, perbc, sorted = true);
