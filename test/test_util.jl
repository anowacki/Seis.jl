using Test
using Seis

@testset "Utility" begin
    # @chain macro
    let b = 0, delta = 1, t = Trace(b, delta, [0,1])
        local f, g
        Seis.@chain f(t::Trace) = t.b + 1
        @test f(t) == (t |> f()) == b + 1
        Seis.@chain function g(t::Trace, x)
            t.b + x
        end
        @test g(t, 2) == (t |> g(2)) == b + 2
    end

    @testset "Nearest sample" begin
        let t = Trace(-2, 1, rand(10))
            @test nearest_sample(t, 3.1) == 6
            @test nearest_sample(t, -5) == 1
            @test nearest_sample(t, -5, inside=true) == nothing
        end
        let t = Trace(0.01, 0.01, rand(100)), time = rand()
            @test nearest_sample(t, time) == round(Int, 100*time)
        end
    end
end
