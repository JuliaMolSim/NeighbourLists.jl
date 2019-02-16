

using NeighbourLists
using JuLIP


# at = set_pbc!(bulk(:Si, cubic=true) * 3, true)
# length(at)
# cutoff = 3.1 * rnn(:Si)
# nlist = CellList(positions(at), cutoff, cell(at), pbc(at))
#
# f(r, R) = (exp(-8*(r-1)) - 2 * exp( - 4 * (r - 1))) * 1 / (1+r)
#
# out = zeros(length(at))
# mapreduce_sym!(out, f, pairs(nlist))
#
# @time mapreduce_sym!(out, f, pairs(nlist));
# @time NeighbourLists.tmapreduce_sym!(out, f, pairs(nlist));



macro symm(N, ex)
   if N isa Symbol
      N = eval(N)
   end
   @assert ex.head == :for
   @assert length(ex.args) == 2
   ex_for = ex.args[1]
   ex_body = ex.args[2]
   # iteration symbol
   i = ex_for.args[1]
   # lower and upper bound
   a0 = ex_for.args[2].args[1]
   a1 = ex_for.args[2].args[2]
   # create the for-loop without body
   loopstr = "for $(i)1 = ($a0):(($a1)-$(N-1))"
   for n = 2:N
      loopstr *= ", $i$n = $i$(n-1)+1:(($a1)-$(N-n))"
   end
   loopstr *= "\n $i = CartesianIndex($(i)1"
   for n = 2:N
      loopstr *= ", $i$n"
   end
   loopstr *= ") \n end"
   loopex = Meta.parse(loopstr)
   append!(loopex.args[2].args, ex_body.args)
   # return the expression
   esc(quote
      $loopex
   end)
end



@generated function t(b0, b1, Nval::Val{N}) where N
   quote
      @symm $N for j = b0:b1
         println(j)
      end
   end
end

t(3, 8, Val(3))

M = 4
@symm M for j = 1:6
   println(j)
end


getN(::Val{N}) where N = N::Int
getN(n::Integer) = n
getN(::Type{Val{N}}) where N = N::Int

macro mm(NN)
   @show NN
   @show getN(eval(NN))
end

@mm(4)
@mm(Val(4))

getN(Val{4})

M = 4
@mm 4
@mm M

function m1(::Val{N}) where N
   @show N
   @show typeof(N)
end

m1(Val(4))

@generated function m2(v::T)
   @mm(T)
end




function simplex_lengths(I::SVector{N, TI}, X::AbstractVector{T}) where {
                                          N, TI <: Integer, T <: AbstractFloat}
   N2 = (N*(N-1))÷2
   a = zeros(TI, N2)
   b = zeros(TI, N2)
   s = zeros(T, N2)
   n = 0
   for i = 1:N-1, j = i+1:N
      n += 1
      a[n] = I[i]
      b[n] = I[j]
      s[n] = abs(X[a[n]] - X[b[n]])
   end
   return s, a, b
end

I = SVector(3, 5, 7, 2)
X = rand(10)
simplex_lengths(I, X)


function statalloc(I::SVector{N, TI}) where {N, TI}
   N2 = (N*(N-1)) ÷ 2
   a = zero(MVector{N2, TI})
end

statalloc(I)

function simplex_lengths!(s, a, b, I::SVector{N, TI}, X) where {N, TI <: Integer}
   n = 0
   for i = 1:N-1, j = i+1:N
      n += 1
      a[n] = I[i]
      b[n] = I[j]
      s[n] = norm(X[a[n]] - X[b[n]])
   end
   return SVector(s), SVector(a), SVector(b)
end

function simplex_lengths2!(s, a, b, I::SVector{N, TI}, X) where {N, TI <: Integer}
   n = 0
   for i = 1:N-1, j = i+1:N
      n += 1
      a[n] = I[i]
      b[n] = I[j]
      s[n] = norm(X[a[n]] - X[b[n]])
   end
   return s, a, b
end



function test(slen, M)
   I = SVector(3, 50, 74, 22, 10)
   N = length(I)
   N2 = (N*(N-1)) ÷ 2
   X = rand(SVector{3, Float64}, 100)
   s = zero(MVector{N2, Float64})
   a = zero(MVector{N2, Int})
   b = zero(MVector{N2, Int})
   for m = 1:M
      slen(s, a, b, I, X)
   end
end



@btime test(simplex_lengths!, 10_000)
@btime test(simplex_lengths2!, 10_000)
# @btime simplex_lengths2!($s, $a, $b, $I, $X)
