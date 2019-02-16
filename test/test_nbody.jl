using ForwardDiff, StaticArrays, NeighbourLists
using Test, LinearAlgebra, Printf

# uncomment for testing from editor/file
# include("test_aux.jl")

let rcut = 3.0

fn, fn_d = gen_fnbody(rcut)

function naive_n_body(X::Vector{SVec{T}}, f, M, rcut) where {T}
   N = length(X)
   E = 0.0
   cnt = 0
   start = CartesianIndex(ntuple(_->1, M))
   stop = CartesianIndex(ntuple(_->N, M))
   for j in CartesianIndices(start, stop)
      s = zeros(T, (M*(M-1))÷2); n = 0
      for a = 1:M-1, b = a+1:M
         n += 1
         s[n] = norm(X[j[a]] - X[j[b]])
      end
      if 0 < minimum(s) && maximum(s) <= rcut
         # each of these terms occur factorial(M) times, i.e. for
         # every permutation of the indices!
         E += f(s) / factorial(M)
      end
   end
   return E
end

vecs(V::Vector{T}) where {T} = reinterpret(SVector{3,T}, V, (length(V) ÷ 3,))
mat(V::Vector{SVector{3,T}}) where {T} = reinterpret(T, V, (3, length(V)))

# NOT NEEDED - WE DO FINITE-DIFFERENCE TESTS INSTEAD!
# function grad_naive_n_body(X, f, M, rcut)
#    x = mat(X)[:]
#    g = ForwardDiff.gradient( y -> naive_n_body(vecs(y), f, M, rcut),  x)
#    return vecs(g)
# end


println("--------------------------------------")
println("    Testing NBodyIterator")
println("--------------------------------------")
println("Check that the energy is consistent with a naive implementation")
MM = [2,2,3,3,3,4,4,5]  # body orders
println("   N     Nat    =>   |Emr-Enaive|")
for M in MM
   # create a not-too-large copper cell
   X, C, _ = rand_config(2)
   nat = length(X)

   # assemble energy via neighbourlist and map-reduce
   Emr = n_body(X, fn, M, rcut, C)
   # assemble energy naively
   Enaive = naive_n_body(X, fn, M, rcut)

   println("   $M      $nat    =>   $(abs(Emr - Enaive))")
   @test Emr ≈ Enaive
end


println("--------------------------------------")
println("Finite-difference tests")
MM = [2,3,4,5]  # body orders
for M in MM
   println(" [M = $M]")
   X, C, pbc = rand_config(2)
   if M > 3
      pbc = [false, false, false]
   end
   @show pbc
   nat = length(X)
   dE = grad_n_body(X, fn_d, M, rcut, C, pbc)
   dE = mat(dE)[:]
   E = n_body(X, fn, M, rcut, C, pbc)
   @printf("    h    | error \n")
   @printf("---------|----------- \n")
   x = mat(X)[:]
   errors = []
   for p = 2:11
      h = 0.1^p
      dEh = copy(dE)
      for n = 1:length(dE)
         x[n] += h
         dEh[n] = (n_body(vecs(x), fn, M, rcut, C, pbc) - E) / h
         x[n] -= h
      end
      push!(errors, vecnorm(dE - dEh, Inf))
      @printf(" %1.1e | %4.2e  \n", h, errors[end])
   end
   @test minimum(errors) <= 1e-3 * maximum(errors)
end


end # let block
