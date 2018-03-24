using ForwardDiff, StaticArrays, NeighbourLists

include("test_aux.jl")


# MODEL N-Body function
rcut = 2.1
fnbody = let rcut = rcut
   rs -> sqrt(sum(exp.(0.5-rs))) .* prod( (rs/rcut-1.0).^2 .* (rs .< rcut) )
end
fnbody_d(rs) = ForwardDiff.gradient(fnbody, rs)

function naive_M_body(X, f, M, rcut)
   N = length(X)
   E = 0.0
   start = CartesianIndex(ntuple(_->1, M))
   stop = CartesianIndex(ntuple(_->N, M))
   for j in CartesianRange(start, stop)
      s = zeros((M*(M-1))÷2); n = 0
      for a = 1:M-1, b = a+1:M
         n += 1
         s[n] = norm(X[j[a]] - X[j[b]])
      end
      if maximum(s) < rcut
         E += f(s) / M
      end
   end
   return E
end

# # check that fnbody works as expected
# rs = @SVector rand(5)
# fnbody(rs)
# fnbody_d(rs)

println("--------------------------------------")
MM = [3,3,3,4,4] # ,5]  # body orders
for M in MM
   # create a not-too-large copper cell
   X, C, _ = rand_config(2)
   nat = length(X)
   @show M, nat
   # assemble energy via neighbourlist and map-reduce
   nlist = PairList(X, rcut, C, (false, false, false), sorted = true)
   Emr = NeighbourLists.mapreduce_sym!(fnbody, zeros(nat),
                                 NeighbourLists.nbodies(M, nlist)) |> sum

   # assemble energy naively
   Enaive = naive_M_body(X, fnbody, M, rcut)

   @show Emr, Enaive
end
