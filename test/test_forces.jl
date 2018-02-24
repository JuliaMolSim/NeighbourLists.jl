
using NeighbourList
using JuLIP, JuLIP.Potentials, StaticArrays


# ==================================
#     Lennard-Jones Test
# ==================================

function lj_iter(nlist, V)
   Es = zeros(nsites(nlist))
   for (i, _1, r, _2) in pairs(nlist)
      Es[i] += V(r)
   end
   dE = zeros(typeof(nlist.R[1]), nsites(nlist))
   for (i, j, r, R) in pairs(nlist)
      dV = (@D V(r)) * (R/r)
      dE[j] += dV
      dE[i] -= dV
   end
   return Es, dE
end


function lj_mapreduce{T}(nlist::CellList{T}, V)
   Es = mapreduce_sym!(zeros(T, nsites(nlist)),
                      (r,R) -> V(r), pairs(nlist) )
   dE = mapreduce_antisym!(zeros(JVec{T}, nsites(nlist)),
                      (r,R) -> ((@D V(r))/r) * R, pairs(nlist) )
end


println("----------------------------------------")
println("Lennard-Jones Test")
println("----------------------------------------")

r0 = rnn("Fe")
cutoff =  r0 * 2.7
lj = LennardJones(r0, 1.0) * C1Shift(cutoff)

for L in (5, 10, 20)
   at = bulk("Fe", cubic=true) * L
   println("Bulk Fe, Nat = $(length(at))")
   println("CellList construction:")
   @time nlist = CellList(positions(at), cutoff, cell(at), pbc(at))
   @time nlist = CellList(positions(at), cutoff, cell(at), pbc(at))
   println("Energy + Force assembly Iterator:")
   @time lj_iter(nlist, lj)
   @time lj_iter(nlist, lj)
   println("Energy + Force assembly MapReduce:")
   @time lj_mapreduce(nlist, lj)
   @time lj_mapreduce(nlist, lj)
end



# ==================================
#     EAM Test
# ==================================


function eam_iter{T}(nlist::CellList{T}, V)
   Es = zeros(T, nsites(nlist))
   for (i, _1, r, R) in NeighbourList.sites(nlist)
      Es[i] += V(r, R)
   end
   dE = zeros(JVec{T}, nsites(nlist))
   for (i, j, r, R) in NeighbourList.sites(nlist)
      dV = @D V(r, R)
      dE[j] += dV
      dE[i] -= sum(dV)
   end
   return Es, dE
end

function eam_mapreduce{T}(nlist::CellList{T}, V)
   Es = map!(V, zeros(T, nsites(nlist)), NeighbourList.sites(nlist))
   dE = map_cfd!((r, R) -> (@D V(r,R)), zeros(JVec{T}, nsites(nlist)),
                NeighbourList.sites(nlist))
end

println("----------------------------------------")
println("Analytic EAM Test")
println("----------------------------------------")
r0 = rnn("Fe")
cutoff = r0 * 2.7
ϕ = (@analytic r -> 1/r) * C1Shift(cutoff)
ρ = (@analytic r -> exp(-r)) * C1Shift(cutoff)
eam = EAM(ϕ, ρ, @analytic t -> sqrt(t))

for L in (5, 10, 20)
   at = bulk("Fe", cubic=true) * L
   println("Bulk Fe, Nat = $(length(at))")
   println("CellList construction:")
   @time nlist = CellList(positions(at), cutoff, cell(at), pbc(at))
   @time nlist = CellList(positions(at), cutoff, cell(at), pbc(at))
   println("Energy + Force assembly Iterator:")
   @time eam_iter(nlist, eam)
   @time eam_iter(nlist, eam)
   println("Energy + Force assembly MapReduce:")
   @time eam_mapreduce(nlist, eam)
   @time eam_mapreduce(nlist, eam)
end
