
using NeighbourLists: SMat, SVec

function rand_config(N)
   C = SMat( diagm(2.0 + 0.2 * rand(3)) * N )
   X = [ C' * rand(SVec)   for i = 1:ceil(Int, abs(det(C))) ÷ 4 + 2 ]
   pbc = SVec(rand(Bool, 3))
   return X, C, pbc
end
