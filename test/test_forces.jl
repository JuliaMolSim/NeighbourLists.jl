
using NeighbourList
using JuLIP

# ==================================
#     Lennard-Jones Test
# ==================================

function V_lj(r)
   s6 = (1/r)^6
   return s6 * (s6 - 2.0)
end

function dV_lj(r)
   s = 1/r
   s6 = s^6
   return - 12 * s6 * s * (s6 - 1.0)
end

function lj_iter(nlist, rnn)
   Es = zeros(nsites(nlist))
   n = 0
   for (i, _1, r, _2) in pairs(nlist)
      Es[i] += V_lj(r/rnn)
   end
   dE = zeros(typeof(nlist.R[1]), nsites(nlist))
   for (i, j, r, R) in pairs(nlist)
      dV = (dV_lj(r/rnn)/rnn/r) * R
      dE[j] += dV
      dE[i] -= dV
   end
   return Es, dE
end


function lj_mapreduce(nlist, rnn)
   Es = mapreduce( (r,R) -> V_lj(r/rnn), pairs(nlist) )
   dE = mapreduce_d( (r,R) -> (dV_lj(r/rnn)/rnn/r) * R, pairs(nlist) )
end


for L in (3, 10, 30)
   cutoff = rnn("Fe") * 2.9
   at = bulk("Fe", cubic=true) * L
   println("Bulk Fe, Nat = $(length(at))")
   println("CellList construction:")
   @time nlist = CellList(positions(at), cutoff, cell(at), pbc(at))
   @time nlist = CellList(positions(at), cutoff, cell(at), pbc(at))
   println("Energy + Force assembly Iterator:")
   @time lj_iter(nlist, rnn("Fe"))
   @time lj_iter(nlist, rnn("Fe"))
   println("Energy + Force assembly MapReduce:")
   @time lj_mapreduce(nlist, rnn("Fe"))
   @time lj_mapreduce(nlist, rnn("Fe"))
end
