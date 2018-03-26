using ForwardDiff, StaticArrays, NeighbourLists
using Base.Test

# uncomment for testing from editor/file
# include("test_aux.jl")

# # MODEL N-Body function
# rcut = 2.1
# fn = let rcut = rcut
#    rs -> sqrt(sum(exp.(0.5-rs))) .* prod( (rs/rcut-1.0).^2 .* (rs .< rcut) )
# end
# fn_d(rs) = ForwardDiff.gradient(fn, rs)

let rcut = 2.1

fn, fn_d = gen_fnbody(rcut)

function naive_M_body{T}(X::Vector{SVec{T}}, f, M, rcut)
   N = length(X)
   E = 0.0
   cnt = 0
   start = CartesianIndex(ntuple(_->1, M))
   stop = CartesianIndex(ntuple(_->N, M))
   for j in CartesianRange(start, stop)
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

vecs{T}(V::Vector{T}) = reinterpret(SVector{3,T}, V, (length(V) ÷ 3,))
mat{T}(V::Vector{SVector{3,T}}) = reinterpret(T, V, (3, length(V)))


function grad_naive_M_body(X, f, M, rcut)
   x = mat(X)[:]
   g = ForwardDiff.gradient( y -> naive_M_body(vecs(y), f, M, rcut),  x)
   return vecs(g)
end

function M_body(X, f, M, rcut, C)
   nlist = PairList(X, rcut, C, (false, false, false), sorted = true)
   return NeighbourLists.mapreduce_sym!(f, zeros(length(X)),
                                 NeighbourLists.nbodies(M, nlist)) |> sum
end

function grad_M_body(X, df, M, rcut, C)
   nlist = PairList(X, rcut, C, (false, false, false), sorted = true)
   return NeighbourLists.mapreduce_sym_d!(df, zeros(SVec{Float64}, length(X)),
                                 NeighbourLists.nbodies(M, nlist))
end



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
   Emr = M_body(X, fn, M, rcut, C)
   # assemble energy naively
   Enaive = naive_M_body(X, fn, M, rcut)

   println("   $M      $nat    =>   $(abs(Emr - Enaive))")
   @test Emr ≈ Enaive
end


println("--------------------------------------")
println("Finite-difference tests")
MM = [2,3,4,5]  # body orders
for M in MM
   println(" [M = $M]")
   X, C, _ = rand_config(2)
   nat = length(X)
   dE = grad_M_body(X, fn_d, M, rcut, C)
   dE = mat(dE)[:]
   E = M_body(X, fn, M, rcut, C)
   @printf("    h    | error \n")
   @printf("---------|----------- \n")
   x = mat(X)[:]
   errors = []
   for p = 2:11
      h = 0.1^p
      dEh = copy(dE)
      for n = 1:length(dE)
         x[n] += h
         dEh[n] = (M_body(vecs(x), fn, M, rcut, C) - E) / h
         x[n] -= h
      end
      push!(errors, vecnorm(dE - dEh, Inf))
      @printf(" %1.1e | %4.2e  \n", h, errors[end])
   end
   @test minimum(errors) <= 1e-3 * maximum(errors)
end



end # let block
