# Trace operations
using Test
using Dates: Second, now
using Statistics
using Seis
import Seis.SAC

"Test that two traces are the same after removing any filename"
function compare_remove_filename!(t1, t2; rtol=√eps(eltype(t1.t)))
    t1.meta.file = missing
    t2.meta.file = missing
    t1.b == t2.b && t1.delta == t2.delta && isapprox(t1.t, t2.t, rtol=rtol)
end

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
            time_now = now()
            t.evt.time = time_now
            @test cut(t, time_now, time_now + Second(1)) == cut(t, 0, 1)
            add_pick!(t, 1, "Test pick")
            @test cut(t, "Test pick", 0, "Test pick", 0.5) == cut(t, 1, 1.5)
            @test cut(t, "Test pick", 0, 0.5) == cut(t, 1, 1.5)
            # Empty traces
            @test_throws ArgumentError cut(t, endtime(t)+1, endtime(t)+2)
            @test_throws ArgumentError cut(t, starttime(t)-2, starttime(t)-1)
            @test nsamples(cut(t, endtime(t)+1, endtime(t)+2, allowempty=true)) == 0
            t′ = cut(t, -2, -1, allowempty=true)
            @test starttime(t′) == -2
            @test endtime(t′) == -2 - t.delta
        end
    end

    @testset "Decimate" begin
        @test compare_copy_modify_func(decimate, decimate!, 2)
        @test compare_copy_modify_func(decimate, decimate!, 2, antialias=false)
        let n = rand(10:1000), t = Trace(1, 0.5, rand(n))
            # With antialiasing
            @test_throws ArgumentError decimate(t, 0)
            @test decimate(t, 5, antialias=true) == decimate(t, 5)
            @test nsamples(decimate(t, 5)) == length((1:n)[1:5:end])
            @test decimate(t, 3).delta == t.delta*3
            t′ = deepcopy(t)
            decimate!(t′, 4)
            @test t′ == decimate(t, 4)
            # Without antialiasing
            @test_throws ArgumentError decimate(t, 0, antialias=false)
            @test decimate(t, 1, antialias=false) == t
            @test nsamples(decimate(t, 5, antialias=false)) == length((1:n)[1:5:end])
            @test decimate(t, 3, antialias=false).delta == t.delta*3
            t′ = deepcopy(t)
            decimate!(t′, 4, antialias=false)
            @test t′ == decimate(t, 4, antialias=false)
        end
        # Test antialiasing filter
        let fs = 5, delta = 1/5, npts = round(Int, 1000/fs), t = Trace(0, delta, npts)
            t.t = [sin(π*tt/2) + sin(3π*tt)/3 for tt in times(t)]
            decimate!(t, 3)
            @test t.t ≈ [sin(π*tt/2) for tt in times(t)] rtol=0.05
        end
    end

    @testset "Differentiation" begin
        let t = sample_data(), t′ = deepcopy(t)
            for i in (0, 1, 4, 6)
                @test_throws ArgumentError differentiate(t, points=i)
            end

            for i in (2, 3, 5)
                local file = joinpath(@__DIR__, "test_data", "operations",
                    "seis_diff_points_$i.sac")
                local reference_data = read_sac(file)
                @test compare_remove_filename!(differentiate(t, points=i), reference_data)
            end
            differentiate!(t′)
            @test t′ == differentiate(t)
        end
    end

    @testset "Flip" begin
        let t = Trace(0, 1, [0, 1, 0])
            local t′
            @test_throws ArgumentError flip(t)
            t.sta.azi = 10
            t.sta.inc = 80
            t′ = deepcopy(t)
            flip!(t′)
            @test flip(t) == t′
            @test t′.sta.azi == 190
            @test t′.sta.inc == 100
            @test trace(t′) == [0, -1, 0]
        end
    end

    @testset "Integration" begin
        let t = sample_data(), t′ = deepcopy(t)
            @test_throws ArgumentError integrate(t, :weird_method)

            for method in (:trapezium, :rectangle)
                local file = joinpath(@__DIR__, "test_data", "operations",
                    "seis_int_$(method).sac")
                local reference_data = read_sac(file)
                @test compare_remove_filename!(integrate(remove_mean(t), method),
                    reference_data)
            end
            integrate!(t′)
            @test t′ == integrate(t)
        end
    end

    @testset "Remove mean" begin
        let t = Trace(0, 0.01, [i%2 for i in 1:1000]), t′ = deepcopy(t)
            atol = eps(eltype(trace(t)))
            @test mean(trace(remove_mean(t))) ≈ 0.0 atol=atol
            remove_mean!(t)
            @test t == remove_mean(t′)
        end
    end

    @testset "Remove trend" begin
        let t = Trace(0, 0.1, 1:100), t′ = deepcopy(t)
            atol = 1000eps(eltype(trace(t)))
            @test all(isapprox.(trace(remove_trend(t)), 0.0, atol=atol))
            remove_trend!(t)
            @test t == remove_trend(t′)
        end
    end

    @testset "Normalisation" begin
        for (f, f!) in zip((normalise, normalize), (normalise!, normalize!))
            let t = Trace(0, 1, rand(5)), val=rand(1:20), t′ = deepcopy(t)
                @test maximum(abs, trace(f(t, val))) ≈ val
                @test maximum(abs, trace(f(t))) ≈ 1
                f!(t)
                @test t == f(t′)
            end
        end
        let t = Trace(0, 1, rand(5)), t′ = deepcopy(t)
            @test normalise(t) == normalize(t)
            normalise!(t)
            normalize!(t′)
            @test t == t′
        end
    end

    @testset "Taper" begin
        let t = Trace(rand(), rand(), rand(100)), t′ = deepcopy(t)
            @test_throws ArgumentError taper(t, 0.6)
            @test_throws ArgumentError taper(t, -0.1)
            @test_throws ArgumentError taper(t, form=:random_symbol)
            @test trace(taper(t))[1] == trace(taper(t))[end] == 0.0
            for form in (:hamming, :hanning, :cosine)
                @test trace(taper(t))[end÷2] == trace(t)[end÷2]
            end
            taper!(t, 0.3)
            @test t == taper(t′, 0.3)
        end
        let t = Trace(0, 1, fill(1.0, 100))
            @test all(isapprox.(trace(taper(t, 0.01, form=:hamming))[3:end-2], 1.0))
            @test all(isapprox.(trace(taper(t, 0.30, form=:cosine))[31:end-30], 1.0))
        end
    end

    @testset "Envelope" begin
        let t = Trace(-rand(), rand(), fill(1.0f0, 1000)), t′ = deepcopy(t)
            @test trace(envelope(t)) ≈ trace(t)
            envelope!(t)
            @test t == envelope(t′)
        end
        let t = Trace(-rand(), rand(), fill(-rand(), 100)), t′ = deepcopy(t)
            @test trace(envelope(t)) ≈ abs.(trace(t))
        end
        @test compare_copy_modify_func(envelope, envelope!)
    end
end
