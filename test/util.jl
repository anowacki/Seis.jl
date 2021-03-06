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

            @testset "is_east" begin
                @test is_east(Station(inc=90, azi=90))
                @test is_east(Station(inc=90, azi=91), tol=1)
                @test !is_east(Station(inc=90, azi=91))
                @test !is_east(Station(inc=90, azi=91), tol=0.1)
                let t = Trace(0, 1, 0)
                    t.sta = Station(inc=90, azi=90)
                    @test is_east(t)
                    @test is_east(t, tol=0.1)
                    t.sta = Station(inc=90, azi=91)
                    @test !is_east(t)
                    @test is_east(t, tol=1)
                    @test !is_east(t, tol=0.1)
                end
            end

            @testset "is_north" begin
                @test is_north(Station(inc=90, azi=0))
                @test is_north(Station(inc=90, azi=1), tol=1)
                @test is_north(Station(inc=90, azi=359), tol=1)
                @test !is_north(Station(inc=90, azi=1))
                @test !is_north(Station(inc=90, azi=1), tol=0.1)
                @test !is_north(Station(inc=90, azi=359), tol=0.1)
                let t = Trace(0, 1, 0)
                    t.sta = Station(inc=90, azi=0)
                    @test is_north(t)
                    @test is_north(t, tol=0.1)
                    t.sta = Station(inc=90, azi=1)
                    @test !is_north(t)
                    @test is_north(t, tol=1)
                    @test !is_north(t, tol=0.1)
                    t.sta = Station(inc=90, azi=359)
                    @test !is_north(t)
                    @test is_north(t, tol=1)
                end
            end

            @testset "origin_time!" begin
                # Adjust picks
                let t = Trace(0, 1, 2)
                    t.evt.time = DateTime(2000)
                    t.picks.A = 10
                    origin_time!(t, DateTime(2000) + Second(10) + Millisecond(123))
                    @test t.evt.time == DateTime(2000, 1, 1, 0, 0, 10, 123)
                    @test startdate(t) == DateTime(2000)
                    @test starttime(t) == t.b ≈ -10.123
                    @test t.picks.A.time ≈ -0.123
                end
                # Don't adjust picks
                let t = Trace(0, 1, 2)
                    t.evt.time = DateTime(2000)
                    t.picks.A = -1
                    origin_time!(t, DateTime(1999, 12, 31, 23, 59, 59), picks=false)
                    @test t.evt.time == DateTime(1999, 12, 31, 23, 59, 59)
                    @test startdate(t) == DateTime(2000)
                    @test starttime(t) == 1
                    @test t.picks.A.time == -1
                end
                # No event time set
                let b = rand(), t = Trace(b, 1, 2), p = rand()
                    t.picks.A = p, "A"
                    origin_time!(t, DateTime(3000))
                    @test t.evt.time == DateTime(3000)
                    @test starttime(t) == b
                    @test startdate(t) == DateTime(3000) + Millisecond(round(Int, 1000*b))
                    @test t.picks.A.time == p
                    @test t.picks.A.name == "A"
                end
            end

            @testset "origin_time" begin
                let t = Trace(0, 1, 2), time = DateTime(3000, 1, 2, 3, 4, 5, 678),
                        time′ = time + Second(rand(-100:100))
                    t.evt.time = time
                    t.picks.A = (1, "A")
                    @test origin_time(t, time′) == origin_time!(deepcopy(t), time′)
                    @test origin_time(t, time′, picks=false) ==
                        origin_time!(deepcopy(t), time′, picks=false)
                end
            end
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
        let t = sample_data()
            @test nearest_sample(t, t.evt.time) === nothing
            @test nearest_sample(t, t.evt.time, inside=false) == 1
        end
    end

    @testset "nsamples" begin
        let n = 100, t = Trace(-1, 0.5, n)
            t.evt.time = DateTime(3000)
            @test nsamples(t, starttime(t), endtime(t)) == nsamples(t) == n
            @test nsamples(t, -1, -0.51) == 1
            @test nsamples(t, 1, 2) == 3
            @test nsamples(t, 1, 2.5) == 4
            @test nsamples(t, -10, 50) == n
            @test nsamples(t, 2, 1) == 0
            @test nsamples(t, -6, -1.1) == 0
            @test nsamples(t, 50, 100) == 0
            @test nsamples(t, -6, -1) == 1
            @test nsamples(t, startdate(t), enddate(t)) == n
            @test nsamples(t, DateTime(2999), DateTime(2999, 12, 31, 23, 59, 59)) == 1
            @test nsamples(t, DateTime(2000), DateTime(4000)) == n
            @test nsamples(t, DateTime(3000), DateTime(3000) + Millisecond(4500)) == 10
            @test nsamples(t, enddate(t), enddate(t) - Second(1)) == 0
            @test nsamples(t, typemin(DateTime), typemax(DateTime)) == n
            @test nsamples(t, typemin(DateTime), DateTime(3000)) == 3
            @test nsamples(t, DateTime(3000), typemax(DateTime)) == 98
            @test nsamples(t, -Inf, Inf) == n
            @test nsamples(t, -Inf, 0) == 3
            @test nsamples(t, 0, Inf) == 98
            @test nsamples(t, Inf, 1) == 0
            @test nsamples(t, 0, -Inf) == 0
        end
    end
end
