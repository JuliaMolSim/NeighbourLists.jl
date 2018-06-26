
using JuLIP, ASE, StaticArrays, NeighbourLists
using Base.Test

X_Ti = vecs([0.0 5.19374 2.59687 3.8953 1.29843 6.49217 7.7906 12.9843 10.3875 11.6859 9.08904 14.2828; 0.0 0.918131 1.83626 -1.11022e-16 0.918131 1.83626 -2.22045e-16 0.918131 1.83626 0.0 0.918131 1.83626; 0.0 0.0 0.0 2.24895 2.24895 2.24895 0.0 0.0 0.0 2.24895 2.24895 2.24895])
C_Ti = (@SMatrix [15.5812 2.47895 0.0; 0.0 2.75439 0.0; 0.0 0.0 4.49791])

# this is now completely meaningless, we need to switch this test to compare
# against ASE / matscipy
function test_nlist_julip(at, cutoff)
   nlist = PairList(positions(at), cutoff, cell(at), pbc(at))
   py_nlist = ASE.matscipy_nlist(at, cutoff)
   return ((nlist.i == py_nlist.i) &&
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

println("JuLIP Configuration tests:")
for (i, (descr, at, cutoff)) in enumerate(test_configs)
  print("TEST $i: $descr => ")
  println(@test test_nlist_julip(at, cutoff))
end
