using Test
using Seis

"Make sure that copying and in-place versions respectively do not and do modify the trace"
function compare_copy_modify_func(f, f!, args...; kwargs...)
    t = Trace(100rand(), rand(), rand(100))
    t′ = deepcopy(t)
    f(t, args...; kwargs...) == f!(t′, args...; kwargs...) && t != t′
end

@testset "Seis" begin
    include("filtering.jl")
    include("io.jl")
    include("operations.jl")
    include("rotation.jl")
    include("sample_data.jl")
    include("synth.jl")
    include("traveltimes.jl")
    include("types.jl")
    include("util.jl")
end
