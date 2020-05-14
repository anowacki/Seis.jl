using Test, Dates
using Seis

@testset "Utility" begin
    @testset "Angle difference" begin
        @test abs(Seis.angle_difference(0, π, false)) ≈ π atol=eps(Float64)
        @test Seis.angle_difference(0, 160, true) == 160
        @test Seis.angle_difference(0, 10) == 10
        @test Seis.angle_difference(10, 0) == -10
        @test Seis.angle_difference(359, -1) == 0
        @test Seis.angle_difference(359, 0) == 1
    end

    @testset "'Get/Setters'" begin
        let b = rand(), delta = rand(), n = rand(1:1000), v = rand(n),
                t = Trace(b, delta, v), t′ = deepcopy(t)
            @test nsamples(t) == n
            @test starttime(t) == b
            @test endtime(t) == b + (n-1)*delta
            @test trace(t) == v
            @test times(t) ≈ b:delta:(b+(n-1)*delta)
            t.sta.inc = 0
            @test is_vertical(t)
            @test is_vertical(t.sta)
            @test !is_horizontal(t)
            @test !is_horizontal(t.sta)
            t.sta.inc = 90
            @test !is_vertical(t)
            @test !is_vertical(t.sta)
            @test is_horizontal(t)
            @test is_horizontal(t.sta)
            t.sta.inc = 91
            @test !is_vertical(t)
            @test !is_vertical(t.sta)
            @test !is_horizontal(t)
            @test !is_horizontal(t.sta)
            t.sta.inc = t′.sta.inc = 90
            t.sta.azi = 0
            t′.sta.azi = 90
            @test Seis.traces_are_orthogonal(t, t′) && Seis.traces_are_orthogonal(t′, t)
            t.sta.azi = 360rand()
            t′.sta.azi = t.sta.azi - 90rand(1:2:100)
            @test Seis.traces_are_orthogonal(t, t′) && Seis.traces_are_orthogonal(t′, t)
        end
    end

    @testset "Dates and times" begin
        # dates
        let b = rand(), delta = 0.001, n = 340, tr = rand(n), t = Trace(b, delta, tr)
            @test_throws ErrorException dates(t) # No origin time defined
            t.evt.time = DateTime(3004, 2, 29, 08, 47, 21, 400)
            @test length(dates(t)) == 340
            @test dates(t)[1] == t.evt.time + Dates.Millisecond(round(Int, t.b*1e3))
            @test Dates.value(dates(t)[2] - dates(t)[1]) ≈ t.delta*1e3 # in ms
        end
        let b = 0, delta = 1, n = 61, t = Trace(b, delta, rand(n))
            t.evt.time = DateTime(1999, 12, 31, 23, 59, 0)
            @test dates(t)[end] == DateTime(2000)
            t.delta = 0.25e-3
            @test_throws ErrorException dates(t) # delta < 1 ms
            @test_throws ErrorException startdate(t)
            @test_throws ErrorException enddate(t)
            t.delta = 1
            t.b = 1
            @test dates(t)[1] == t.evt.time + Dates.Second(1)
            @test dates(t)[1] == startdate(t)
            @test dates(t)[end] == DateTime(2000, 1, 1, 0, 0, 1)
            @test dates(t)[end] == enddate(t)
        end
    end

    @testset "@chain" begin
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

    @testset "Nearest sample" begin
        let t = Trace(-2, 1, rand(10))
            @test nearest_sample(t, 3.1) == 6
            @test nearest_sample(t, -5, inside=false) == 1
            @test nearest_sample(t, -5) == nothing
        end
        let t = Trace(0, 0.01, rand(101)), time = rand()
            @test nearest_sample(t, time) == round(Int, 100*time) + 1
        end
        let t = Trace(0, 1, rand(61))
            t.evt.time = DateTime(2000, 1, 1)
            @test nearest_sample(t, DateTime(2000, 1, 1, 0, 0, 30)) == 31
            @test nearest_sample(t, DateTime(3000), inside=false) == nsamples(t)
            @test nearest_sample(t, DateTime(1000)) == nothing
        end
    end
end
