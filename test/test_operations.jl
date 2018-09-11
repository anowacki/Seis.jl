# Trace operations
using Compat.Test
import Dates: Second, now
using Seis
import SAC

@testset "Operations" begin
    @testset "Cut" begin
        let t = Trace(0, 0.01, rand(200)), t′ = deepcopy(t), b = rand(), e = b+rand()
            @test_throws ArgumentError cut(t, e, b)
            @test_throws ArgumentError cut(t, missing, e)
            @test cut(t, t.b, endtime(t)) == t
            @test cut(t, t.b-1, endtime(t)+1, warn=false) == t
            @test_logs (:warn,
                "Beginning cut time -1.0 is before start of trace.  Setting to 0.0.")
                (:warn,
                "End cut time 2.99 is after end of trace.  Setting to 1.99.")
                cut(t, t.b-1, endtime(t)+1)
            @test_logs cut(t, t.b-1, endtime(t)+1, warn=false)
            @test cut!(t′, b, e) == cut(t, b, e)
            @test t′ == cut(t, b, e)
            @test t′.b ≈ b atol=t.delta/2
            @test endtime(t′) ≈ e atol=t.delta/2
            t.evt.time = now()
            @test cut(t, now(), now() + Second(1)) == cut(t, 0, 1)
            add_pick!(t, 1, "Test pick")
            @test cut(t, "Test pick", 0, "Test pick", 0.5) == cut(t, 1, 1.5)
            @test cut(t, "Test pick", 0, 0.5) == cut(t, 1, 1.5)
        end
    end
end