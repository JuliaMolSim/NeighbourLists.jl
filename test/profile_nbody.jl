using NeighbourLists
using JuLIP
using BenchmarkTools

include("test_aux.jl")

print("N-body energy and forces benchmark: ")
println("# Threads = ", Base.Threads.nthreads())

at = bulk(:Cu, cubic=true) * 5
println("bulk :Cu with nat = $(length(at))")
r0 = rnn(:Cu)
rcut = 2.1 * r0
X = positions(at)
C = cell(at)
f, f_d = gen_fnbody(rcut, r0)
nlist = PairList(X, rcut, C, (false, false, false), sorted = true)

for M in [2, 3, 4, 5]
   println("M = $M")
   print("  energy: ")
   @btime n_body($f, $M, $nlist)
   print("  forces: ")
   @btime grad_n_body($f_d, $M, $nlist)
end


lj = LennardJones(r0, 1.0) * C1Shift(cutoff)
@btime forces($lj, $at)
