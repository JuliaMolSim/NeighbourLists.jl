using ForwardDiff, StaticArrays, NeighbourLists
using Base.Test

# uncomment for testing from editor/file
# include("test_aux.jl")

# MODEL N-Body function
rcut = 2.1
fnbody = let rcut = rcut
   rs -> sqrt(sum(exp.(0.5-rs))) .* prod( (rs/rcut-1.0).^2 .* (rs .< rcut) )
end
fnbody_d(rs) = ForwardDiff.gradient(fnbody, rs)

function naive_M_body(X, f, M, rcut)
   N = length(X)
   E = 0.0
   cnt = 0
   start = CartesianIndex(ntuple(_->1, M))
   stop = CartesianIndex(ntuple(_->N, M))
   for j in CartesianRange(start, stop)
      s = zeros((M*(M-1))÷2); n = 0
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

# # check that fnbody works as expected
# rs = @SVector rand(5)
# fnbody(rs)
# fnbody_d(rs)


println("--------------------------------------")
println("    Testing NBodyIterator")
println("--------------------------------------")
MM = [2,2,3,3,3,4,4,5]  # body orders
println("   N     Nat    =>   |Emr-Enaive|")
for M in MM
   # create a not-too-large copper cell
   X, C, _ = rand_config(2)
   nat = length(X)

   # assemble energy via neighbourlist and map-reduce
   nlist = PairList(X, rcut, C, (false, false, false), sorted = true)
   Emr = NeighbourLists.mapreduce_sym!(fnbody, zeros(nat),
                                 NeighbourLists.nbodies(M, nlist)) |> sum

   # assemble energy naively
   Enaive = naive_M_body(X, fnbody, M, rcut)

   println("   $M      $nat    =>   $(abs(Emr - Enaive))")
   @test abs(Emr - Enaive) < 1e-10
end
