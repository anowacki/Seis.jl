using Test
using Seis

@testset "Seis" begin
    include("filtering.jl")
    include("io.jl")
    include("operations.jl")
    include("sample_data.jl")
    include("synth.jl")
    include("traveltimes.jl")
    include("types.jl")
    include("util.jl")
end
