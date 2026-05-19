using NeighbourLists
using Test

# Include shared test utilities
include("test_utils.jl")

# ---- FLAGS -----

# whether to run performance tests
performance = true

# check whether on CI
isCI = haskey(ENV, "CI")
notCI = !isCI

# Initialize GPU backend detection (supports CUDA, AMDGPU, Metal)
init_gpu_backend!()

# Legacy compatibility
cuda_available = check_cuda_available()

# ----------------- TESTS -------------------

println("# threads = $(Base.Threads.nthreads())")
println("# GPU backend = $(gpu_backend())")

@testset "NeighbourLists" begin
   @testset "CellList (Legacy)" begin include("test_celllist.jl") end
   @testset "SortBased" begin include("test_sortbased.jl") end
   @testset "Unified API" begin include("test_unified_api.jl") end
   @testset "AtomsBase Integration" begin include("test_atoms_base.jl") end
   @testset "GPU" begin include("test_gpu.jl") end
   @testset "Matscipy" begin include("test_matscipy.jl") end
end
