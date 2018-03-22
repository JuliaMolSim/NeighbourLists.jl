
using JuLIP
using Base.Test

if hasjulip

   # this is now completely meaningless, we need to switch this test to compare
   # against ASE / matscipy
   function test_nlist_julip(at, cutoff)
      cl = PairList(positions(at), cutoff, cell(at), pbc(at))
      nlist = neighbourlist(at, cutoff)
      return (cl.i == nlist.i) && (cl.j == nlist.j)
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
   ]

   println("JuLIP Configuration tests:")
   for (i, (descr, at, cutoff)) in enumerate(test_configs)
      println("TEST $i: $descr")
      @test test_nlist_julip(at, cutoff)
   end

end
