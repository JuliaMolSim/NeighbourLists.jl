# This file contains some code that is chared across several tests

using StaticArrays, Test, ForwardDiff
using NeighbourLists: SMat, SVec
using Test, StaticArrays, ForwardDiff

# ------ generate random configurations -------


function rand_config(N)
   C = SMat( diagm(0 => 2.0 .+ 0.2 * rand(3)) * N )
   X = [ C' * rand(SVec)   for i = 1:ceil(Int, abs(det(C))) ÷ 4 + 2 ]
   pbc = SVec(rand(Bool, 3))
   return X, C, pbc
end


# --------- MANY BODY CODE THAT IS SHARED ACROSS TESTS ------------


sqrt1p(x) = expm1( 0.5 * log1p( x ) )   # sqrt(1 + x) - 1
sqrt1p_d(x) = 0.5 / sqrt(1+x)

function fnbody(r, r0, rcut)
   E = sum( exp.(1.0 .- r./r0) )
   F = prod( (r./rcut - 1.0) .* (r .< rcut) )
   return sqrt1p(E) * F^2
end

function fnbody_d(r, r0, rcut)
   E = sum( exp.(1.0 .- r./r0) )
   F = prod( (r./rcut - 1.0) .* (r .< rcut) )
   ∇f = (- sqrt1p_d(E) * F^2 / r0) .* (exp.(1.0 .- r./r0))
   ∇f += (sqrt1p(E) * 2.0 / rcut * F^2) .* (r./rcut .- 1.0 - 1e-14).^(-1)
   return ∇f
end

# Generate a MODEL N-Body function
function gen_fnbody(rcut, r0=1.0)
   return r->fnbody(r, r0, rcut), r -> fnbody_d(r, r0, rcut)
end

println("Checking that `fnbody` is correct...")
fnbody_ad(rs, r0, rcut) = ForwardDiff.gradient( t -> fnbody(t, r0, rcut), rs )
for n = 1:5
   r = rand(SVector{3, Float64})
   @test fnbody_d(r, 1.0, 2.0) ≈ fnbody_ad(r, 1.0, 2.0)
end

# global assembly of n-body energies and forces

n_body(X, f, M, rcut, C, pbc = (false, false, false),
            nlist = PairList(X, rcut, C, pbc, sorted = true)) =
   n_body(f, M, nlist)

n_body(f, M, nlist::PairList) =
   NeighbourLists.maptosites!(f, zeros(nsites(nlist)),
                                 NeighbourLists.nbodies(M, nlist)) |> sum

grad_n_body(X, df, M, rcut, C, pbc = (false, false, false),
            nlist = PairList(X, rcut, C, pbc, sorted = true)) =
   grad_n_body(df, M, nlist)

grad_n_body(df, M, nlist::PairList) =
   NeighbourLists.maptosites_d!(df, zeros(SVec{Float64}, nsites(nlist)),
                                  NeighbourLists.nbodies(M, nlist))
