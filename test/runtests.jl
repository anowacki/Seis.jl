using Compat.Test
using Seis

@testset "Seis" begin
    include("test_types.jl")
    include("test_io.jl")
    include("test_traveltimes.jl")
    include("test_util.jl")
end