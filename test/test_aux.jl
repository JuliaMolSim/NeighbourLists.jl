
using NeighbourLists: SMat, SVec

function rand_config(N)
   C = SMat( diagm(2.0 + 0.2 * rand(3)) * N )
   X = [ C' * rand(SVec)   for i = 1:ceil(Int, abs(det(C))) ÷ 4 + 2 ]
   pbc = SVec(rand(Bool, 3))
   return X, C, pbc
end


# --------- MANY BODY CODE THAT IS SHARED ACROSS TESTS ------------


sqrt1p(x) = expm1( 0.5 * log1p( x ) )

function fnbody(r, r0, rcut)
   E = 0.0
   F = 1.0
   for n = 1:length(r)
      E += exp(1.0 - r[n]/r0)
      F *= (r[n]/rcut - 1.0) * (r[n] < rcut)
   end
   return sqrt1p(E) * F^2
end

# fnbody(rs, r0, rcut) = sqrt(sum(exp.(1.0-rs/r0))) .* prod( (rs/rcut-1.0).^2 .* (rs .< rcut) )

# function fnbody_d(r, r0, rcut)
#    E = 0.0
#    for n = 1:length(r)
#       E += exp(1.0 - r[n]/r0) * (r[n]/rcut - 1.0)^2
#    end
#    rtE = sqrt(1+E)
#    return (0.5 / sqrt(1+E)) * (
#               (-1/r0) * exp.(1.0 .- r./r0) .* (r/rcut .- 1.0).^2
#             + exp.(1.0 .- r./r0) .* (2.0/rcut) .* (r./rcut .- 1.0) )
# end

fnbody_d(rs, r0, rcut) = ForwardDiff.gradient( t -> fnbody(t, r0, rcut), rs )


# Generate a MODEL N-Body function
function gen_fnbody(rcut, r0=1.0)
   # return r->fnbody(r, r0, rcut), r -> fnbody_d(r, r0, rcut)
   return r->fnbody(r, r0, rcut), r -> fnbody_d(r, r0, rcut)
end

#
# fnbody(rs, r0, rcut) =
#    sqrt(sum(exp.(0.5-rs))) .* prod( (rs/rcut-1.0).^2 .* (rs .< rcut) )
#
# fnbody_d(r, r0, rcut) = ForwardDiff.gradient(r->fnbody(r,r0,rcut), r)
#
#
# # println("Checking that `fnbody` is correct...")
# # for n = 1:5
# #    r = rand(SVector{3, Float64})
# #    @test fnbody_d(r, 1.0, 2.0) ≈ fnbody_ad(r, 1.0, 2.0)
# # end
#
#
# # global assembly of n-body energies and forces
#
# n_body(X, f, M, rcut, C,
#             nlist = PairList(X, rcut, C, (false, false, false), sorted = true)) =
#       n_body(f, M, nlist)
#
# n_body(f, M, nlist::PairList) =
#    NeighbourLists.mapreduce_sym!(f, zeros(nsites(nlist)),
#                                  NeighbourLists.nbodies(M, nlist)) |> sum
#
# grad_n_body(X, df, M, rcut, C,
#             nlist = PairList(X, rcut, C, (false, false, false), sorted = true)) =
#    grad_n_body(df, M, nlist)
#
# grad_n_body(df, M, nlist::PairList) =
#    NeighbourLists.mapreduce_sym_d!(df, zeros(SVec{Float64}, nsites(nlist)),
#                                   NeighbourLists.nbodies(M, nlist))
