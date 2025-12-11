using NeighbourLists
using Test

# ---- FLAGS -----

# whether to run performance tests
performance = true

# check whether on CI
isCI = haskey(ENV, "CI")
notCI = !isCI

# Check for CUDA availability (for GPU tests)
cuda_available = false
try
   using CUDA
   global cuda_available = CUDA.functional()
catch
   # CUDA not available
end

# ----------------- TESTS -------------------

println("# threads = $(Base.Threads.nthreads())")
println("# CUDA available = $(cuda_available)")

@testset "NeighbourLists" begin
   @testset "Aux" begin include("test_aux.jl") end
   @testset "CellList (Legacy)" begin include("test_celllist.jl") end
   @testset "SortBased" begin include("test_sortbased.jl") end
   include("test_atoms_base.jl")
   @testset "GPU" begin include("test_gpu.jl") end
   @testset "Matscipy" begin include("test_matscipy.jl") end
end
