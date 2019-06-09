
using JuLIP, ASE, StaticArrays, NeighbourLists
using Test

X_Ti = vecs([0.0 5.19374 2.59687 3.8953 1.29843 6.49217 7.7906 12.9843 10.3875 11.6859 9.08904 14.2828; 0.0 0.918131 1.83626 -1.11022e-16 0.918131 1.83626 -2.22045e-16 0.918131 1.83626 0.0 0.918131 1.83626; 0.0 0.0 0.0 2.24895 2.24895 2.24895 0.0 0.0 0.0 2.24895 2.24895 2.24895])
C_Ti = (@SMatrix [15.5812 2.47895 0.0; 0.0 2.75439 0.0; 0.0 0.0 4.49791])



# -------------- MatSciPy NeighbourList Patch -------------
using PyCall
import NeighbourLists
matscipy_neighbours = pyimport("matscipy.neighbours")
function asenlist(at::Atoms, rcut)
   pyat = ASEAtoms(at).po
   return matscipy_neighbours[:neighbour_list]("ijdD", pyat, rcut)
end

function matscipy_nlist(at::Atoms{T}, rcut::T; recompute=false, kwargs...) where T <: AbstractFloat
   i, j, r, R = asenlist(at, rcut)
   i = copy(i) .+ 1
   j = copy(j) .+ 1
   r = copy(r)
   R = collect(vecs(copy(R')))
   first = NeighbourLists.get_first(i, length(at))
   NeighbourLists.sort_neigs!(j, r, R, first)
   return NeighbourLists.PairList(positions(at), rcut, i, j, r, R, first)
end
# --------------------------------------------------------

pynlist(at, cutoff) = matscipy_nlist(at, cutoff)
jnlist(at, cutoff) = PairList(positions(at), cutoff, cell(at), pbc(at))

function test_nlist_julip(at, cutoff)
   nlist = jnlist(at, cutoff)
   py_nlist = pynlist(at, cutoff)
   return ( (nlist.i == py_nlist.i) &&
            (nlist.j == py_nlist.j) &&
            (nlist.r ≈ py_nlist.r) &&
            (nlist.R ≈ py_nlist.R) )
end


test_configs = [
   #
   ( "si, cubic, cluster, short",
    set_pbc!(bulk(:Si, cubic=true) * 3, false),
    1.1 * rnn(:Si) ),
   #
   ( "si, cubic, cluster, med",
    set_pbc!(bulk(:Si, cubic=true) * 3, false),
    2.1 * rnn(:Si) ),
   #
   ( "si, non-cubic cell, cluster, med",
    set_pbc!(bulk(:Si) * 5, false),
    2.1 * rnn(:Si) ),
   #
   ( "si, non-cubic cell, pbc",
    set_pbc!(bulk(:Si) * 5, false),
    2.1 * rnn(:Si) ),
   #
   ( "si, non-cubic cell, mixed bc",
    set_pbc!(bulk(:Si) * 5, false),
    2.1 * rnn(:Si) ),
    #
   ("Ti, non-symmetric elongated cell, pbc",
    set_pbc!( set_cell!( Atoms(:Ti, X_Ti), C_Ti ), true ),
    1.45 * rnn(:Ti) ),
   #
   ("Ti (hcp?) canonical cell, pbc",
    set_pbc!( bulk(:Ti) * 4, true ),
    2.3 * rnn(:Ti) ),
   ]

# ----------- A FEW MORE COMPLEX TESTS THAT FAILED AT SOME POINT ---------------

# [1] a left-handed cell orientation
#     the test that failed during Cas' experiments

X = [ 0.00000000e+00  0.00000000e+00  0.00000000e+00
      1.92333044e+00  6.63816518e-17 -1.36000000e+00
      1.92333044e+00  1.92333044e+00 -2.72000000e+00
      3.84666089e+00  1.92333044e+00 -4.08000000e+00 ]'
C = [ 3.84666089   0.           0.
      0.           3.84666089   0.
      0.           0.          -5.44 ]'

# X = [ 0.00000000e+00  0.00000000e+00  0.00000000e+00
#       1.92333044e+00  6.63816518e-17  1.36000000e+00
#       1.92333044e+00  1.92333044e+00  2.72000000e+00
#       3.84666089e+00  1.92333044e+00  4.08000000e+00 ]'
# C = [ 3.84666089   0.           0.
#       0.           3.84666089   0.
#       0.           0.           5.44 ]'
at = Atoms(:Si, collect(vecs(X)))
set_cell!(at, C)
set_pbc!(at, (true,true,true))
atlge = at * (1,1,10)
rcut = 2.3*rnn(:Si)
push!(test_configs, ("Si left-oriented", at, rcut))
push!(test_configs, ("Si left-oriented, large", atlge, rcut))

test_nlist_julip(at, rcut)
test_nlist_julip(atlge, rcut)

# # [2] vacancy in bulk Si
#
# using PyCall
# at = bulk(:Si, cubic=true)
# @pyimport ase.lattice.cubic as cubic
# at2py = cubic.Diamond(symbol = "Si", latticeconstant = 5.43)
# at2py[:get_positions]()
# @test at2py[:get_cell]() ≈ cell(at1)
# @test mat(positions(at1))' ≈ at2py[:get_positions]()[:,[3,2,1]]
#
# at = bulk(:Si, cubic=true)
# at1 = at * 3
# X = mat(positions(at1))
# at1 = deleteat!(at1, 1)
# at2 = deleteat!(set_positions!(at * 3, X[[3,2,1],:]), 1)
# rcut = 2.3 * rnn(:Si)
# test_nlist_julip(at1, rcut)
# test_nlist_julip(at2, rcut)

# [3] Two Ti configuration that seems to be causing problems
C1 = @SMatrix [5.71757 -1.81834e-15 9.74255e-41; -2.85879 4.95156 4.93924e-25; 4.56368e-40 9.05692e-25 9.05629]
X1 = vecs([0.00533847 2.85879 -1.42939 1.42939 0.0 2.85879 -1.42939 1.42939 -1.43e-6 2.85878 -1.42939 1.42939 -1.43e-6 2.85878 -1.42939 1.42939;
          -0.0 -0.0 2.47578 2.47578 0.0 -0.0 2.47578 2.47578 1.65052 1.65052 4.1263 4.1263 1.65052 1.65052 4.1263 4.1263;
          0.00845581 0.0 0.0 0.0 4.52815 4.52815 4.52815 4.52815 2.26407 2.26407 2.26407 2.26407 6.79222 6.79222 6.79222 6.79222]) |> collect
C2 = @SMatrix [5.71757 0.0 0.0; -2.85879 4.95156 0.0; 0.0 0.0 9.05629]
X2 = vecs([0.00534021 2.85879 -1.42939 1.42939 0.0 2.85879 -1.42939 1.42939 0.0 2.85879 -1.4294 1.4294 0.0 2.85879 -1.4294 1.4294;
           0.0 0.0 2.47578 2.47578 0.0 0.0 2.47578 2.47578 1.65052 1.65052 4.1263 4.1263 1.65052 1.65052 4.1263 4.1263;
           0.00845858 0.0 0.0 0.0 4.52815 4.52815 4.52815 4.52815 2.26407 2.26407 2.26407 2.26407 6.79222 6.79222 6.79222 6.79222])  |> collect
at1 = set_cell!(Atoms(:Ti, X1), C1)
at2 = set_cell!(Atoms(:Ti, X2), C2)
rcut = 2.5 * rnn(:Ti)
test_nlist_julip(at1, rcut)
test_nlist_julip(at2, rcut)

# --------------- ACTUALLY RUNNING THE TESTS ------------------

println("JuLIP Configuration tests:")
for (i, (descr, at, cutoff)) in enumerate(test_configs)
  print("TEST $i: $descr => ")
  println(@test test_nlist_julip(at, cutoff))
end
