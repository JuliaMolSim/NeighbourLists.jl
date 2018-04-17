
using NeighbourLists
using JuLIP, StaticArrays


# ==================================
#     Lennard-Jones Test
# ==================================

function lj_iter(nlist, V)
   Es = zeros(nsites(nlist))
   for (i, _1, r, _2) in pairs(nlist)
      Es[i] += V(r)
   end
   dE = zeros(eltype(nlist.R), nsites(nlist))
   for (i, j, r, R) in pairs(nlist)
      dV = (@D V(r)) * (R/r)
      dE[j] += dV
      dE[i] -= dV
   end
   return Es, dE
end


function lj_mapreduce{T}(nlist::PairList{T}, V)
   Es = maptosites!( (r,R) -> V(r), zeros(T, nsites(nlist)),
                        pairs(nlist) )
   dE = maptosites_d!((r,R) -> (@D V(r)), zeros(JVec{T}, nsites(nlist)),
                         pairs(nlist) )
end

function lj_nbody{T}(nlist::PairList{T}, V)
   Es = maptosites!(r -> V(r[1]),  zeros(T, nsites(nlist)),
                      nbodies(2, nlist) )
   dE = maptosites_d!(r -> (@D V(r[1])),  zeros(JVec{T}, nsites(nlist)),
                      nbodies(2, nlist) )
end


println("----------------------------------------")
println("Lennard-Jones Test")
println("----------------------------------------")

r0 = rnn(:Fe)
cutoff =  r0 * 2.7
lj = LennardJones(r0, 1.0) * C1Shift(cutoff)

for L in (2, 5, 10, 20)
   at = bulk(:Fe, cubic=true) * L
   println("Bulk Fe, Nat = $(length(at))")
   println("PairList construction:")
   @time nlist = PairList(positions(at), cutoff, cell(at), pbc(at))
   @time nlist = PairList(positions(at), cutoff, cell(at), pbc(at))
   println("Energy + Force assembly Iterator:")
   @time lj_iter(nlist, lj)
   @time lj_iter(nlist, lj)
   println("Energy + Force assembly MapReduce:")
   @time lj_mapreduce(nlist, lj)
   @time lj_mapreduce(nlist, lj)
   println("Energy + Force assembly nbody MapReduce:")
   @time lj_nbody(nlist, lj)
   @time lj_nbody(nlist, lj)
end



# # ==================================
# #     EAM Test
# # ==================================
#
#
# function eam_iter{T}(nlist::PairList{T}, V)
#    Es = zeros(T, nsites(nlist))
#    for (i, _1, r, R) in NeighbourLists.sites(nlist)
#       Es[i] += V(r, R)
#    end
#    dE = zeros(JVec{T}, nsites(nlist))
#    for (i, j, r, R) in NeighbourLists.sites(nlist)
#       dV = @D V(r, R)
#       dE[j] += dV
#       dE[i] -= sum(dV)
#    end
#    return Es, dE
# end
#
# function eam_iter!{T}(nlist::PairList{T}, V)
#    Es = zeros(T, nsites(nlist))
#    for (i, _1, r, R) in NeighbourLists.sites(nlist)
#       Es[i] += V(r, R)
#    end
#    dE = zeros(JVec{T}, nsites(nlist))
#    ndv = maximum(length(j) for (i, j, r, R) in NeighbourLists.sites(nlist))
#    dV = zeros(JVec{T}, ndv)
#    for (i, j, r, R) in NeighbourLists.sites(nlist)
#       fill!(dV, zero(JVec{T}))
#       JuLIP.Potentials.evaluate_d!(dV, V, r, R)
#       for a = 1:length(j)
#          dE[j[a]] += dV[a]
#          dE[i] -= dV[a]
#       end
#    end
#    return Es, dE
# end
#
# function eam_mapreduce{T}(nlist::PairList{T}, V)
#    Es = map!(V, zeros(T, nsites(nlist)), NeighbourLists.sites(nlist))
#    dE = map_cfd!((r, R) -> (@D V(r,R)), zeros(JVec{T}, nsites(nlist)),
#                 NeighbourLists.sites(nlist))
# end
#
# println("----------------------------------------")
# println("Analytic EAM Test")
# println("----------------------------------------")
# r0 = rnn("Fe")
# cutoff = r0 * 2.7
# ϕ = (@analytic r -> 1/r) * C1Shift(cutoff)
# ρ = (@analytic r -> exp(-r/3)) * C1Shift(cutoff)
# eam = EAM(ϕ, ρ, @analytic t -> sqrt(1+t))
#
# for L in (5, 10, 20)
#    at = bulk("Fe", cubic=true) * L
#    println("Bulk Fe, Nat = $(length(at))")
#    println("PairList construction:")
#    @time nlist = PairList(positions(at), cutoff, cell(at), pbc(at))
#    @time nlist = PairList(positions(at), cutoff, cell(at), pbc(at))
#    println("Energy + Force assembly Iterator:")
#    @time eam_iter(nlist, eam)
#    @time eam_iter(nlist, eam)
#    println("Energy + Force assembly MapReduce:")
#    @time eam_mapreduce(nlist, eam)
#    @time eam_mapreduce(nlist, eam)
#    println("Energy + Force assembly In-place Iterator:")
#    @time eam_iter!(nlist, eam)
#    @time eam_iter!(nlist, eam)
# end
