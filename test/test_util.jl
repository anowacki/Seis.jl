using Compat.Test
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
end
