
using NeighbourLists, StaticArrays
using JuLIP

include("test_configs.jl")

# _, at, cutoff = test_configs[2]

at = bulk("Si", cubic=true) * 10
set_pbc!(at, (true, false, true))
cutoff = 2.1 * rnn("Si")

i, j, r, R = nn_list(positions(at), cutoff, SMatrix{3,3}(cell(at)),
                        SVec(pbc(at)...), one(Int))

clist = CellList(positions(at), cutoff, cell(at), pbc(at), sorted = true)

i == clist.i
j == clist.j

@time nn_list(positions(at), cutoff, SMatrix{3,3}(cell(at)),
                        SVec(pbc(at)...), one(Int))
@time nn_list(positions(at), cutoff, SMatrix{3,3}(cell(at)),
                        SVec(pbc(at)...), one(Int))
@time nn_list(positions(at), cutoff, SMatrix{3,3}(cell(at)),
                        SVec(pbc(at)...), one(Int))


@time  CellList(positions(at), cutoff, cell(at), pbc(at), sorted = true)
@time  CellList(positions(at), cutoff, cell(at), pbc(at), sorted = true)
@time  CellList(positions(at), cutoff, cell(at), pbc(at), sorted = true)


#
# SUMMARY:
#   These tests should demonstrate that NearestNeighbors.jl is faster
#   for small configurations but not for large ones. And this is even
#   before general cell shapes are taken into account.
#


using BenchmarkTools

@btime nn_list(positions(at), cutoff, SMatrix{3,3}(cell(at)),
                        SVec(pbc(at)...), one(Int))


@btime  CellList(positions(at), cutoff, cell(at), pbc(at), sorted = true)
