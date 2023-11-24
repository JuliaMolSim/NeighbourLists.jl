using NeighbourLists
using JuLIP
using BenchmarkTools

##


# using PyCall
# using ASE
# using ProfileView

# include("../test/test_aux.jl")

# @pyimport matscipy.neighbours as matscipy_neighbours
# matscipy_neighbours = pyimport("matscipy.neighbours")

# function matscipy_nlist(at, cutoff)
#    pycall(matscipy_neighbours["neighbour_list"],
#           NTuple{5, PyArray}, "ijdDS", at.po, cutoff)
# end

println("# Threads = ", Base.Threads.nthreads())

##

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

   # println("Julia Nlist")
   @btime PairList($X, $cutoff, $C, $perbc)
   # println("Matscipy Nlist")
   # @btime matscipy_nlist($(ASEAtoms(at)), $cutoff)
   println("------------------------------------------")
end

##

# at = bulk(:Si, cubic=true, pbc = true ) * 10
# rcut = 2.5 * rnn(:Si)
# C = JMat(cell(at))
# X = positions(at)
# perbc = JVec(pbc(at))

# # NeighbourLists._pairlist_(X, C, perbc, rcut, Int, false)
# clist = NeighbourLists._celllist_(X, C, perbc, rcut, Int)

# @btime NeighbourLists._celllist_($X, $C, $perbc, $rcut, Int)
# @btime NeighbourLists._pairlist_($clist)

# # @btime NeighbourLists._pairlist_($X, $C, $perbc, $rcut, Int, false)
# #    clist = _celllist_(X, cell, pbc, cutoff, TI)
# #    i, j, S = _pairlist_(clist)

# @code_warntype NeighbourLists._pairlist_(clist)

##

# @profview let clist = clist 
#    for _ = 1:100 
#       NeighbourLists._pairlist_(clist)
#    end
# end

##

# at = bulk(:Si, cubic=true, pbc = true ) * 10
# rcut = 2.5 * rnn(:Si)

# @profview let C = JMat(cell(at)), X = positions(at), perbc = JVec(pbc(at)), 
#               rcut = rcut 
#    for _ = 1:30 
#       PairList(X, rcut, C, perbc)
#    end
# end