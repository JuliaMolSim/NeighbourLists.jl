

using NeighbourLists
using JuLIP


at = set_pbc!(bulk(:Si, cubic=true) * 3, true)
length(at)
cutoff = 3.1 * rnn(:Si)
nlist = CellList(positions(at), cutoff, cell(at), pbc(at))

f(r, R) = (exp(-8*(r-1)) - 2 * exp( - 4 * (r - 1))) * 1 / (1+r)

out = zeros(length(at))
mapreduce_sym!(out, f, pairs(nlist))

@time mapreduce_sym!(out, f, pairs(nlist));
@time NeighbourLists.tmapreduce_sym!(out, f, pairs(nlist));
