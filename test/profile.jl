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

# println("Neighbourlist assembly benchmark")
# for L in [4, 10, 30]
#    print("L = $L")
#    # si, non-cubic cell, mixed bc
#    at = bulk(:Si, cubic=true) * L
#    println(", N = $(length(at))")
#    set_pbc!(at, (true, false, true))
#    C = JMat(cell(at))
#    X = positions(at)
#    perbc = JVec(pbc(at))
#    cutoff = 3.5 * rnn(:Si)
#
#    println("Julia Nlist")
#    @btime PairList($X, $cutoff, $C, $perbc, sorted = false)
#    println("Matscipy Nlist")
#    @btime matscipy_nlist($(ASEAtoms(at)), $cutoff)
#    println("------------------------------------------")
# end

print("N-body energy and forces benchmark: ")
at = bulk(:Cu, cubic=true) * 5
println("bulk :Cu with nat = $(length(at))")
r0 = rnn(:Cu)
rcut = 2.2 * r0
X = positions(at)
C = cell(at)
f, f_d = gen_fnbody(rcut, r0)
nlist = PairList(X, rcut, C, (false, false, false), sorted = true)

for M in [2, 3, 4]
   println("M = $M")
   print("  energy: ")
   @btime M_body($f, $M, $nlist)
   print("  forces: ")
   @btime grad_M_body($f_d, $M, $nlist)
end


# # @profile M_body(f, 3, nlist)
# # using ProfileView
# # Base.Profile.print()
# #
# # it = NeighbourLists.nbodies(3, nlist)
# # out = zeros(length(X))
# # @code_warntype NeighbourLists.mapreduce_sym!(f, out, it)
# # @code_warntype NeighbourLists.mt_split_interlaced(100)
# # @code_warntype NeighbourLists.mr_sym_inner!(f, out, it, 1:10)
#
# using StaticArrays
# r = rand(SVector{3, Float64})
# @code_warntype fnbody_d(r, 1.0, 2.0)
#
#
# fnbody_ad(r, r0, rcut) = ForwardDiff.gradient(r->fnbody(r,r0,rcut), r)
# fnbody_d(r, 0.78, 2.12) - fnbody_ad(r, 0.78, 2.12)
