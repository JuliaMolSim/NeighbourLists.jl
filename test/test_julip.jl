
using JuLIP, ASE, StaticArrays, NeighbourLists
using Base.Test

X_Ti = vecs([0.0 5.19374 2.59687 3.8953 1.29843 6.49217 7.7906 12.9843 10.3875 11.6859 9.08904 14.2828; 0.0 0.918131 1.83626 -1.11022e-16 0.918131 1.83626 -2.22045e-16 0.918131 1.83626 0.0 0.918131 1.83626; 0.0 0.0 0.0 2.24895 2.24895 2.24895 0.0 0.0 0.0 2.24895 2.24895 2.24895])
C_Ti = (@SMatrix [15.5812 2.47895 0.0; 0.0 2.75439 0.0; 0.0 0.0 4.49791])

pynlist(at, cutoff) = ASE.matscipy_nlist(at, cutoff)
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
at = Atoms(:Si, vecs(X))
set_cell!(at, C)
set_pbc!(at, (true,true,true))
atlge = at * (1,1,10)
rcut = 2.3*rnn(:Si)
push!(test_configs, ("Si left-oriented", at, rcut))
push!(test_configs, ("Si left-oriented, large", atlge, rcut))

# [2] vacancy in bulk Si

# using PyCall
# at = bulk(:Si, cubic=true)
# @pyimport ase.lattice.cubic as cubic
# at2py = cubic.Diamond(symbol = "Si", latticeconstant = 5.43)
# at2py[:get_positions]()
# @test at2py[:get_cell]() ≈ cell(at1)
# @test mat(positions(at1))' ≈ at2py[:get_positions]()[:,[3,2,1]]

# at = bulk(:Si, cubic=true)
# at1 = at * 3
# X = mat(positions(at1))
# at1 = deleteat!(at1, 1)
# at2 = deleteat!(set_positions!(at * 3, X[[3,2,1],:]), 1)
# rcut = 2.3 * rnn(:Si)
# test_nlist_julip(at1, rcut)
# test_nlist_julip(at2, rcut)


# --------------- ACTUALLY RUNNING THE TESTS ------------------

println("JuLIP Configuration tests:")
for (i, (descr, at, cutoff)) in enumerate(test_configs)
  print("TEST $i: $descr => ")
  println(@test test_nlist_julip(at, cutoff))
end
