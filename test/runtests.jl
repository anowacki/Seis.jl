using Test
using Seis

@testset "Seis" begin
    include("test_filtering.jl")
    include("test_io.jl")
    include("test_operations.jl")
    include("test_sample_data.jl")
    include("test_synth.jl")
    include("test_traveltimes.jl")
    include("test_types.jl")
    include("test_util.jl")
end
