using Test, Dates
using NanoDates: NanoDates, NanoDate
using Seis
using LinearAlgebra: ×, ⋅
import Rotations

using .TestHelpers

@testset "Utility" begin
    @testset "Angle difference" begin
        @test abs(Seis.angle_difference(0, π, false)) ≈ π atol=eps(Float64)
        @test Seis.angle_difference(0, 160, true) == 160
        @test Seis.angle_difference(0, 10) == 10
        @test Seis.angle_difference(10, 0) == -10
        @test Seis.angle_difference(359, -1) == 0
        @test Seis.angle_difference(359, 0) == 1

        @testset "Promotion" begin
            let Ts = (Int, Float32, Float64)
                @testset "$T1" for T1 in Ts
                    α = one(T1)
                    for T2 in Ts
                        β = one(T2)
                        Δ = zero(float(promote_type(T1, T2)))
                        @test Seis.angle_difference(α, β) === Δ
                    end
                end
            end
        end
    end

    @testset "_angle_tol" begin
        @testset "Not values" begin
            @test_throws ArgumentError Seis._angle_tol(1)
        end
        @testset "Not ints" begin
            @test_throws MethodError Seis._angle_tol(Int)
        end
        @testset "$T" for T in (Float16, Float32, Float64)
            t = Trace{T,Vector{T},Seis.Geographic{T}}(0, 1, 0)
            sta = Station{T}()
            tol = T == Float64 ? 1000*√eps(T) :
                  T == Float16 ? 10*√eps(T) :
                  √eps(T)
            @test Seis._angle_tol(T) === tol
            @test Seis._angle_tol(t) === tol
            @test Seis._angle_tol(sta) === tol
        end
        @testset "Multiple arguments" begin
            @test Seis._angle_tol(Float16) == Seis._angle_tol(Float32, Float64, Float16)
            @test Seis._angle_tol(Trace{Float32}(0, 1, 0), Trace{Float64}(0, 1, 0)) ==
                Seis._angle_tol(Float32)
        end
    end

    @testset "_direction_vector" begin
        for (azi, inc, correct) in ((0, 90, [0, 1, 0]), (45, 90, [√2/2, √2/2, 0]),
                (0, 0, [0, 0, 1]))
            @test Seis._direction_vector(azi, inc) ≈ correct
            s = Station(azi=azi, inc=inc)
            t = Trace(0, 1, 0)
            t.sta = s
            @test Seis._direction_vector(t) == Seis._direction_vector(s) ==
                Seis._direction_vector(azi, inc)
        end
    end

    @testset "_direction_to_azimuth_incidence" begin
        atol = 1e-6
        azi, inc = Seis._direction_to_azimuth_incidence([1, 0, 0])
        @test azi ≈ 90 atol=atol
        @test inc ≈ 90 atol=atol
        azi, inc = Seis._direction_to_azimuth_incidence([0, -1, 0])
        @test azi ≈ 180 atol=atol
        @test inc ≈ 90 atol=atol
        azi, inc = Seis._direction_to_azimuth_incidence([0, 0, -10])
        @test inc ≈ 180 atol=atol
    end

    @testset "_directions_are_orthogonal" begin
        @testset "$T" for T in (Float16, Float32, Float64)
            tol = Seis._angle_tol(T)
            v, v′ = vector_and_deviated_vector(90)

            @test !Seis._directions_are_orthogonal(v, v, 0)

            @test Seis._directions_are_orthogonal(v, v′, tol)
            @test Seis._directions_are_orthogonal(v, v′, 0.1)
            @test Seis._directions_are_orthogonal(v, v′, tol)

            w, w′ = vector_and_deviated_vector(91)
            @test !Seis._directions_are_orthogonal(w, w′, tol)
            @test Seis._directions_are_orthogonal(w, w′, 2)
        end
    end

    @testset "_directions_are_parallel" begin
        @testset "$T" for T in (Float16, Float32, Float64)
            tol = Seis._angle_tol(T)
            v, v′ = vector_and_deviated_vector(1)

            @test Seis._directions_are_parallel(v, v, 0)
            @test Seis._directions_are_parallel(v, v, tol)
            @test !Seis._directions_are_parallel(v, v′, 0)
            @test !Seis._directions_are_parallel(v, v′, tol)
            @test Seis._directions_are_parallel(v, v′, 2)
            @test !Seis._directions_are_parallel(v, -v, tol)
        end
    end

    @testset "_directions_are_antiparallel" begin
        @testset "$T" for T in (Float16, Float32, Float64)
            tol = Seis._angle_tol(T)
            v, v′ = vector_and_deviated_vector(179)

            @test Seis._directions_are_antiparallel(v, -v, 0)
            @test Seis._directions_are_antiparallel(v, -v, tol)
            @test !Seis._directions_are_antiparallel(v, v′, tol)
            @test !Seis._directions_are_antiparallel(v, v′, 0)
            @test !Seis._directions_are_antiparallel(v, v′, tol)
            @test Seis._directions_are_antiparallel(v, v′, 2)
            @test !Seis._directions_are_antiparallel(v, -v′, tol)
        end
    end

    @testset "_u_dot_v_and_theta" begin
        @testset "Known" begin
            udotv, θ = Seis._u_dot_v_and_theta([1, 0, 0], [0, -1, 0])
            @test udotv ≈ 0
            @test θ ≈ 90

            udotv, θ = Seis._u_dot_v_and_theta([1, 0, 0], [-1, 0, 0])
            @test udotv ≈ -1
            @test θ ≈ 180
        end

        @testset "Random" begin
            θ = 10
            u, v = vector_and_deviated_vector(θ)
            udotv = sum(x -> x[1]*x[2], zip(u, v))
            udotv′, θ′ = Seis._u_dot_v_and_theta(u, v)
            @test udotv ≈ udotv′
            @test θ ≈ θ′
        end
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

            @testset "are_orthogonal" begin
                @testset "Horizontals" begin
                    @test are_orthogonal(t.sta, t′.sta)
                    @test are_orthogonal(t′.sta, t.sta)
                    @test are_orthogonal(t, t′)
                    @test are_orthogonal(t′, t)
                end

                @testset "Random orientations" begin
                    @testset "$T" for T in (Float16, Float32, Float64)
                        x, y, z = random_basis_traces(T)
                        # x component rotated about vector between y and z by 5°
                        x′ = deepcopy(x)
                        x⃗ = Seis._direction_vector(x)
                        y⃗ = Seis._direction_vector(y)
                        z⃗ = Seis._direction_vector(z)
                        k⃗ = normalize(y⃗ + z⃗)
                        x⃗′ = Seis._rotate_by_vector(x⃗, k⃗, deg2rad(5))
                        azi, inc = Seis._direction_to_azimuth_incidence(x⃗′)
                        x′.sta.azi, x′.sta.inc = azi, inc

                        @testset "Pairs" begin
                            @test are_orthogonal(x, y)
                            @test are_orthogonal(y, x)
                            @test are_orthogonal(y, z)
                            @test are_orthogonal(z, y)
                            @test are_orthogonal(z, x)
                            @test are_orthogonal(x, z)

                            # Different to default tolerance
                            @test !are_orthogonal(x′, y)
                            @test !are_orthogonal(x′, z)
                            @test !are_orthogonal(y, x′)
                            @test !are_orthogonal(z, x′)

                            # Same within 10°
                            @test are_orthogonal(x′, y, tol=10)
                            @test are_orthogonal(x′, z, tol=10)
                            @test are_orthogonal(y, x′, tol=10)
                            @test are_orthogonal(z, x′, tol=10)
                        end

                        @testset "$([traces...].sta.cha)" for traces in trace_permutations(x, y, z)
                            @test are_orthogonal(traces...)
                        end
                    end
                end
            end

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
                    @test startdate(t) == NanoDate(3000) + Nanosecond(round(Int, 1_000_000_000*b))
                    @test t.picks.A.time == p
                    @test t.picks.A.name == "A"
                end
            end

            @testset "origin_time" begin
                let t = Trace(0, 1, 2), time = DateTime(3000, 1, 2, 3, 4, 5, 678),
                        time′ = time + Second(rand(-100:100))
                    t.evt.time = time
                    t.picks.A = (1, "A")
                    @test origin_time(t) == time
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
            @test dates(t)[1] == t.evt.time + Dates.Nanosecond(round(Int, t.b*1e9))
            @test Dates.value(dates(t)[2] - dates(t)[1]) ≈ t.delta*1e9 # in ns
        end
        let b = 0, delta = 1, n = 61, t = Trace(b, delta, rand(n))
            t.evt.time = NanoDate(1999, 12, 31, 23, 59, 0)
            @test dates(t)[end] == DateTime(2000)
            t.delta = 0.25e-9
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

        @testset "datetimes" begin
            n = rand(10:100)
            b = 100*(rand() - 0.5)
            delta = rand() + 0.01
            t = Trace(b, delta, n)
            @test_throws ErrorException datetimes(t)
            origin_time!(t, NanoDate("2000-01-01T01:23:45.012345678"))
            dts = datetimes(t)
            @test dts isa Vector{DateTime}
            @test length(dts) == nsamples(t)
            @test all(i -> dts[i] <= dts[i+1], eachindex(dts)[begin:(end-1)])
        end

        @testset "dates v datetimes" begin
            n = rand(10:100)
            b = 100rand()
            delta = rand() + 0.01
            t = Trace(b, delta, n)
            origin_time!(t, NanoDate("2000-01-01T01:23:45.012345678"))
            @test length(dates(t)) == length(datetimes(t)) == nsamples(t)
            @test first(datetimes(t)) == DateTime(startdate(t))
            @test first(datetimes(t)) == DateTime(
                origin_time(t) + Dates.Nanosecond(floor(Int64, b*1000000000))
            )
            @test last(datetimes(t)) == DateTime(enddate(t))
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

    @testset "_quantise" begin
        @test Seis._quantise(-5, 2, 0.4) == -5.6
        @test Seis._quantise(100, 1.34, 0.22) ≈ 99.38 atol=eps(100.0)
        @test Seis._quantise(-5, 2, -0.1) == -4.1
    end
end
