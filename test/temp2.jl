using ResumableFunctions
using BenchmarkTools

@resumable function g(N)::Vector{Int}
   for n = 1:N
      @yield rand(-10:10, 20)
   end
   Int[]
end

f_resumable(N) = sum( sum(x) for x in g(N) )

function f_loop(N)
   s = 0
   for n = 1:N
      s += sum(rand(-10:10, 20))
   end
   return s
end

f_iterator(N) = sum( sum(rand(-10:10, 20)) for x in 1:N )

N = 100_000
@btime f_resumable(N)
@btime f_loop(N)
@btime f_iterator(N)
