using Test
using Seis

include("TestHelpers.jl")
using .TestHelpers

@testset "Seis" begin
    include("filtering.jl")
    include("io.jl")
    include("operations.jl")
    include("geometry.jl")
    include("rotation.jl")
    include("sample_data.jl")
    include("synth.jl")
    include("traveltimes.jl")
    include("types.jl")
    include("util.jl")

    @testset "IO" begin
        include("io/sac.jl")
        include("io/miniseed.jl")
    end

    include("plot.jl")
end
