using NeighbourLists
using JuLIP
using BenchmarkTools
using PyCall
using ASE
# using ProfileView

include("test_aux.jl")

# @pyimport matscipy.neighbours as matscipy_neighbours
matscipy_neighbours = pyimport("matscipy.neighbours")

function matscipy_nlist(at, cutoff)
   pycall(matscipy_neighbours["neighbour_list"],
          NTuple{5, PyArray}, "ijdDS", at.po, cutoff)
end

println("# Threads = ", Base.Threads.nthreads())

println("Neighbourlist assembly benchmark")
for L in [2, 4, 10, 30]
   print("L = $L")
   # si, non-cubic cell, mixed bc
   at = bulk(:Si, cubic=true) * L
   println(", N = $(length(at))")
   set_pbc!(at, (true, true, true))
   C = JMat(cell(at))
   X = positions(at)
   perbc = JVec(pbc(at))
   cutoff = 3.5 * rnn(:Si)

   println("Julia Nlist")
   @btime PairList($X, $cutoff, $C, $perbc, sorted = false)
   println("Matscipy Nlist")
   @btime matscipy_nlist($(ASEAtoms(at)), $cutoff)
   println("------------------------------------------")
end

print("N-body energy and forces benchmark: ")
at = bulk(:Cu, cubic=true) * 10
println("bulk :Cu with nat = $(length(at))")
r0 = rnn(:Cu)
rcut = 2.1 * r0
X = positions(at)
C = cell(at)
f, f_d = gen_fnbody(rcut, r0)
nlist = PairList(X, rcut, C, (false, false, false), sorted = true)

for M in [2, 3] #, 4, 5]
   println("M = $M")
   print("  energy: ")
   @btime n_body($f, $M, $nlist)
   print("  forces: ")
   @btime grad_n_body($f_d, $M, $nlist)
end
