# Trace operations
using Test
using Dates: DateTime, Second, now
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
            @testset "Missing arguments $f" for f in (cut, cut!)
                @test_throws ArgumentError f(t, e, b)
                @test_throws ArgumentError f(t, missing, e)
                @test_throws ArgumentError f(t, b, missing)
                @test_throws ArgumentError f(t, missing, missing)
            end

            @testset "Times" begin
                @test cut(t, t.b, endtime(t)) == t
                @test cut(t, t.b-1, endtime(t)+1, warn=false) == t
                @test_logs(
                    (:warn,
                    "Beginning cut time -1.0 is before start of trace.  Setting to 0.0."),
                    (:warn,
                    "End cut time 2.99 is after end of trace.  Setting to 1.99."),
                    cut(t, t.b-1, endtime(t)+1)
                )
                @test_logs cut(t, t.b-1, endtime(t)+1, warn=false)
                @test cut!(t′, b, e) == cut(t, b, e)
                @test t′ == cut(t, b, e)
                @test t′.b ≈ b atol=t.delta/2
                @test endtime(t′) ≈ e atol=t.delta/2
            end

            @testset "Dates" begin
                time_now = now()
                t.evt.time = time_now
                @test cut(t, time_now, time_now + Second(1)) == cut(t, 0, 1)
                t′ = deepcopy(t)
                @test cut!(t′, time_now, time_now + Second(1)) ==
                    cut(t, time_now, time_now + Second(1))
                @test t′.b == t.b
                @test startdate(t′) == time_now
                @test enddate(t′) == time_now + Second(1)
            end

            @testset "Picks" begin
                add_pick!(t, 1, "Test pick")
                @test cut(t, "Test pick", 0, "Test pick", 0.5) == cut(t, 1, 1.5)
                @test cut(t, "Test pick", 0, 0.5) == cut(t, 1, 1.5)
            end

            @testset "Empty traces" begin
                @test_throws ArgumentError cut(t, endtime(t)+1, endtime(t)+2)
                @test_throws ArgumentError cut(t, starttime(t)-2, starttime(t)-1)
                @test nsamples(cut(t, endtime(t)+1, endtime(t)+2, allowempty=true)) == 0
                t′ = cut(t, -2, -1, allowempty=true)
                @test starttime(t′) == -2
                @test endtime(t′) == -2 - t.delta
            end
        end

        # Fail gracefully when no picks match
        @testset "Picks" begin
            @testset "$f" for f in (cut, cut!)
                types = (Symbol, Regex, String)
                @testset "Type $T1" for T1 in types
                    pick1 = T1("nopick1")
                    let t = Trace(0, 1, 10)
                        @test_throws ArgumentError f(t, pick1, 1, 2)
                        @testset "Type $T2" for T2 in types
                            pick2 = "nopick2"
                            @test_throws ArgumentError f(t, pick1, 1, pick2, 2)
                            t.picks.A = (1, "A")
                            @test_throws ArgumentError f(t, pick1, 1, T2("A"), 2)
                            @test_throws ArgumentError f(t, T2("A"), 1, pick2, 2)
                        end
                        @test_throws ArgumentError f(t, pick1, 1, 2)
                    end
                end
            end
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

    @testset "Resample" begin
        t = Trace(0, 0.5, sin.(0:0.3:6π) .+ 0.2sin.(0:1.2:24π))
        @testset "Calling" begin
            @test_throws ArgumentError resample(t)
            @test_throws ArgumentError resample(t, n=1, delta=t.delta)
            @test_throws ArgumentError resample!(t)
            @test_throws ArgumentError resample!(t, n=1, delta=t.delta)

        end
        @testset "n = $n" for n in (0.5, 1//2, 1, 4, 8//2, 10)
            @test resample(t, n=n) == resample!(deepcopy(t), n=n)
            @test resample(t, n=n).delta ≈ t.delta/n
            # Check that the points line up approximately
            if n > 1
                @test trace(t) ≈ trace(resample(t, n=n))[1:Int(n):end] atol=1e-3
            end
            # Check that nothing is done if not needed
            if n == 1
                @test trace(resample(t, n=n)) == trace(t)
                @test trace(resample!(t, n=n)) === trace(t)
            end
        end
        @testset "delta = $delta" for delta in (1, 0.5, 0.25, 0.125, 0.05)
            @test resample(t, delta=delta) == resample!(deepcopy(t), delta=delta)
            @test resample(t, delta=delta).delta ≈ delta
            # Check that nothing is done if not needed
            if delta == t.delta
                @test trace(resample(t, delta=delta)) == trace(t)
                @test trace(resample!(t, delta=delta)) === trace(t)
            end
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
        @testset "Empty traces" begin
            let t = Trace(0, 1, zeros(5))
                @test trace(normalise(t)) == trace(t)
                @test trace(normalise(t, rand())) == trace(t)
            end
        end
    end

    @testset "Taper" begin
        let t = Trace(rand(), rand(), rand(100)), t′ = deepcopy(t)
            @test_throws ArgumentError taper(t, 0.6)
            @test_throws ArgumentError taper(t, -0.1)
            @test_throws ArgumentError taper(t, form=:random_symbol)
            @test_throws ArgumentError taper(t, time=-1)
            @test_throws ArgumentError taper(t, time=t.delta*100)
            @test trace(taper(t))[1] == trace(taper(t))[end] == 0.0
            for form in (:hamming, :hanning, :cosine)
                @test trace(taper(t))[end÷2] == trace(t)[end÷2]
                @test trace(taper(t, left=false))[1:end÷2] == trace(t)[1:end÷2]
                @test trace(taper(t, right=false))[end÷2:end] == trace(t)[end÷2:end]
                @test trace(taper(t, left=false, right=false)) == trace(t)
                @test taper(t, 0.01) == taper(t, time=0.01*nsamples(t)*t.delta)
            end
            taper!(t, 0.3)
            @test t == taper(t′, 0.3)
        end

        @testset "One-sided" begin
            let t = Trace(0, 1, fill(1.0, 100))
                @test all(isapprox.(trace(taper(t, 0.01, form=:hamming))[3:end-2], 1.0))
                @test all(isapprox.(trace(taper(t, 0.30, form=:cosine))[31:end-30], 1.0))
            end

            @testset "Width <= 1" begin
                let t = Trace(rand(), rand(), rand(100)), t′ = deepcopy(t)
                    @test_throws ArgumentError taper(t, 1.1; left=false)
                    @test_throws ArgumentError taper(t; left=false, time=t.delta*(nsamples(t) + 1))
                    @testset "Right" begin
                        @test trace(taper(t, 0.8; left=false))[1:20] == trace(t)[1:20]
                        @test trace(taper(t, 0.8; left=false))[21] != trace(t)[21]
                        @test trace(taper(t, 0.8; left=false))[end] == 0
                    end
                    @testset "Left" begin
                        @test trace(taper(t, 0.8; right=false))[end-19:end] == trace(t)[end-19:end]
                        @test trace(taper(t, 0.8; right=false))[end-20] != trace(t)[end-20]
                        @test trace(taper(t, 0.8; right=false))[begin] == 0
                    end
                end
            end
        end

        @testset "_taper_core!" begin
            @testset "n too large" begin
                @test_throws BoundsError Seis._taper_core!(rand(3), 4, true, false, :hanning)
            end

            @testset "Matrices" begin
                @testset "$form" for form in (:hamming, :hanning, :cosine)
                    v = rand(100)
                    m = [v v v]
                    n = 20
                    Seis._taper_core!(v, n, true, true, form)
                    Seis._taper_core!(m, n, true, true, form)
                    for col in eachcol(m)
                        @test v == col
                    end
                end
            end
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

    @testset "Merging" begin
        @testset "merge[!]" begin
            @testset "Incompatiable args" begin
                t = Trace(0, 1, [])
                @testset "$f" for f in (merge, merge!)
                    @test_throws ArgumentError f(t, t; gaps=:linear, taper=0.5)
                end
            end

            @testset "Mismatched channel codes" begin
                t1 = Trace(0, 1, [1, 2, 3])
                t2 = Trace(3, 1, [4, 5, 6])
                t1.sta.net = t2.sta.net = "AB"
                t1.sta.sta = t2.sta.sta = "BCD"
                t1.sta.loc = t2.sta.loc = "00"
                t1.sta.cha = t2.sta.cha = "HHZ"
                for f in (:net, :sta, :loc, :cha)
                    t2′ = deepcopy(t2)
                    setproperty!(t2′.sta, f, "")
                    @test_throws ArgumentError merge!(deepcopy(t1), t2′)
                    @test_throws ArgumentError merge(t1, t2′)
                    true_t = Trace(0, 1, [1, 2, 3, 4, 5, 6])
                    true_t.sta = deepcopy(t1.sta)
                    @test merge!(deepcopy(t1), t2′; check=false) == true_t
                    @test merge!(deepcopy(t1), t2′; check=false) == merge(t1, t2′; check=false)
                end
            end

            @testset "Different sampling rates" begin
                @test_throws ArgumentError merge(Trace(0, 1, []), Trace(0, 2, []))
                @test_throws ArgumentError merge!(Trace(0, 1, []), Trace(0, 2, []))
            end

            @testset "sample_tol" begin
                t1 = Trace(10.1, 0.1, [1, 2, 3])
                t2 = Trace(10.52, 0.1, [5, 6, 7])
                @testset "Incorrect value" begin
                    @test_throws ArgumentError merge(t1, t2; sample_tol=-0.1)
                    @test_throws ArgumentError merge(t1, t2; sample_tol=0.6)
                end

                @test_throws ArgumentError merge(t1, t2)
                @test_throws ArgumentError merge(t1, t2; sample_tol=0.01)
                @test merge(t1, t2; sample_tol=0.5) == Trace(10.1, 0.1, [1, 2, 3, 0, 5, 6, 7])
            end

            @testset "Empty traces" begin
                @testset "All empty" begin
                    @test merge(Trace(0, 1, []), Trace(2, 1, []), Trace(3, 1, [])) == Trace(0, 1, [])
                end

                @testset "First empty" begin
                    @test merge(
                        Trace(-1, 1, []),
                        [
                            Trace(0, 1, [1, 2, 3]),
                            Trace(1, 1, [2, 3, 4])
                        ]
                    ) == Trace(0, 1, [1, 2, 3, 4])

                    @test merge(
                        Trace(-1, 1, []),
                        Trace(2, 1, [1, 2, 3])
                    ) == Trace(2, 1, [1, 2, 3])
                end

                @testset "Others empty" begin
                    @test merge(
                        Trace(0, 1, [1, 2, 3]),
                        [
                            Trace(-1, 1, []),
                            Trace(4, 1, [5, 6])
                        ]
                    ) == Trace(0, 1, [1, 2, 3, 0, 5, 6])
                end
            end

            @testset "Overlap options" begin
                @testset "Perfect overlap" begin
                    t1 = Trace(0, 1, [1, 2, 3])
                    t2 = Trace(0, 1, [2, 3, 4])
                    t3 = Trace(0, 1, [6, 7, 8])
                    @testset ":mean" begin
                        @test merge(t1, t2, t3) == Trace(0, 1, [3, 4, 5])
                        @test merge(t1, t2, t3; overlaps=:mean) == Trace(0, 1, [3, 4, 5])
                    end
                    @testset ":first" begin
                        @test merge(t1, t2, t3; overlaps=:first) == t1
                    end
                    @testset ":last" begin
                        @test merge(t1, t2, t3; overlaps=:last) == t3
                    end
                    @testset ":error" begin
                        @test_throws Exception merge(t1, t2, t3; overlaps=:error)
                    end
                    @testset "Invalid option" begin
                        @test_throws ArgumentError merge(t1, t2, t3; overlaps=:not_an_option)
                    end
                end

                @testset "Partial overlap" begin
                    t1 = Trace(-5.5, 0.5, [1, 2, 3, 4])
                    t2 = Trace(-4.5, 0.5, [11, 12, 13, 14])
                    @test merge(t1, t2) == Trace(-5.5, 0.5, [1, 2, 7, 8, 13, 14])
                    @test merge(t1, t2; overlaps=:mean) == Trace(-5.5, 0.5, [1, 2, 7, 8, 13, 14])
                    @test merge(t1, t2; overlaps=:first) == Trace(-5.5, 0.5, [1, 2, 3, 4, 13, 14])
                    @test merge(t1, t2; overlaps=:last) == Trace(-5.5, 0.5, [1, 2, 11, 12, 13, 14])

                    @test merge(t2, t1) == Trace(-5.5, 0.5, [1, 2, 7, 8, 13, 14])
                    @test merge(t2, t1; overlaps=:mean) == Trace(-5.5, 0.5, [1, 2, 7, 8, 13, 14])
                    @test merge(t2, t1; overlaps=:first) == Trace(-5.5, 0.5, [1, 2, 3, 4, 13, 14])
                    @test merge(t2, t1; overlaps=:last) == Trace(-5.5, 0.5, [1, 2, 11, 12, 13, 14])
                end
            end

            @testset "Gap options" begin
                t1 = Trace(0, 1, [0, 1, 2, 3])
                t2 = Trace(-5, 1, [-5, -4, -3, -2])
                t3 = Trace(10, 1, [10, 5, 0])

                @testset "Default :zero" begin
                    @test merge(t1, t2, t3) == merge(t1, t2, t3; gaps=:zero)
                end

                @testset ":zero" begin
                    @test merge(t1, t2, t3) == Trace(-5, 1,
                        [-5, -4, -3, -2, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 10, 5, 0])
                end

                @testset ":error" begin
                    @test_throws Exception merge(t1, t2, t3; gaps=:error)
                end

                @testset ":linear" begin
                    @test merge(t1, t2, t3; gaps=:linear) == Trace(-5, 1,
                        [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 5, 0])
                end

                @testset "value" begin
                    @test merge(t1, t2, t3; gaps=-999) == Trace(-5, 1,
                        [-5, -4, -3, -2, -999, 0, 1, 2, 3, -999, -999, -999, -999, -999, -999, 10, 5 , 0])
                end

                @testset "Invalid option" begin
                    @test_throws ArgumentError merge(t1, t2, t3; gaps=:not_an_option)
                end

                @testset "Tapering" begin
                    @testset "Wrong form" begin
                        @test_throws ArgumentError merge(Trace(0, 1, [1]), Trace(1, 1, [2]), taper_form=:wrong_form)
                    end

                    @testset "Taper" begin
                        t1 = Trace(100, 0.5, randn(10))
                        t2 = Trace(25, 0.5, randn(60))
                        t3 = Trace(0, 0.5, randn(40))

                        @testset "$form" for form in (:hanning, :hamming, :cosine)
                            t_merged = merge(t1, t2, t3; taper=0.2, taper_form=form)

                            @test Trace(0, 0.5, trace(t_merged)[1:40]) ==
                                taper(t3, left=false, form=form, time=0.2*10*0.5)
                            @test Trace(100, 0.5, trace(t_merged)[end-9:end]) ==
                                taper(t1, right=false, form=form, time=4.5)

                            t2′ = deepcopy(t2)
                            taper!(t2′, right=false, form=form, time=0.2*10*0.5)
                            taper!(t2′, left=false, form=form, time=0.2*90*0.5)
                            @test Trace(25, 0.5, trace(t_merged)[51:110]) == t2′
                        end
                    end
                end
            end

            @testset "Sequential, not overlapping" begin
                t1, t2, t3 = t = Trace.([0, 3, 6], 1, Ref([1,2,3]))
                @test merge(t1, t2, t3) == Trace(0, 1, [1, 2, 3, 1, 2, 3, 1, 2, 3])
                @test merge!(t1, t2, t3) == Trace(0, 1, [1, 2, 3, 1, 2, 3, 1, 2, 3])
                @test t1 == Trace(0, 1, [1, 2, 3, 1, 2, 3, 1, 2, 3])
            end

            @testset "Perfect overlap" begin
                @test merge(
                    Trace(5, 0.5, [-2, -3, -4]),
                    Trace(5, 0.5, [2, 3, 4])
                ) == Trace(5, 0.5, [0, 0, 0])
            end

            @testset "Containment" begin
                t1 = Trace(5, 1, [10, 12, 14])
                t2 = Trace(2, 1, [4, 6, 8, 10, 12, 14, 16, 18])
                @test merge(
                    t1, t2
                ) == Trace(2, 1, [4, 6, 8, 10, 12, 14, 16, 18])
            end

            @testset "Overlap" begin
                t1 = Trace(0, 0.01, fill(1.0, 3000))
                t2 = Trace(15, 0.01, fill(3.0, 3000))
                @test merge(t1, t2; overlaps=:mean) == Trace(0, 0.01,
                    [fill(1.0, 1500); fill(2.0, 1500); fill(3.0, 1500)])
            end
        end

        @testset "Absolute merging" begin
            t1 = Trace(0, 1, [0, 1, 2, 3])
            t1.evt.time = DateTime(2000)
            t2 = Trace(0, 1, [3, 4, 5, 6])
            t2.evt.time = DateTime(2000, 1, 1, 0, 0, 1)
            t3 = Trace(0, 1, [-5, -4, -3])
            t3.evt.time = DateTime(1999, 12, 31, 23, 59, 55)

            @testset "relative default" begin
                t_test = Trace(-5, 1, [-5, -4, -3, NaN, NaN, 0, 2, 3, 4, 6])
                t_test.evt.time = t1.evt.time

                @test merge(t1, t2, t3; gaps=NaN, overlaps=:mean) == t_test
            end

            @testset "relative = true" begin
                t_test = Trace(0, 1, [-5, -4, -3, 6])
                t_test.evt.time = t1.evt.time
                @test merge(t1, t2, t3; relative=true, overlaps=:last) == t_test
            end
        end

        @testset "Arrays and tuples" begin
            t1 = Trace(0, 1, [0, 1, 2, 3])
            t1.evt.time = DateTime(2000)
            t2 = Trace(0, 1, [3, 4, 5, 6])
            t2.evt.time = DateTime(2000, 1, 1, 0, 0, 1)
            t3 = Trace(0, 1, [-5, -4, -3])
            t3.evt.time = DateTime(1999, 12, 31, 23, 59, 55)

            @test merge([t1, t2, t3]) == merge(t1, [t2, t3])
            @test merge(t1, t2, t3) == merge(t1, [t2, t3])
            @test merge(t1, (t2, t3)) == merge(t1, [t2, t3])
        end

        @testset "Metadata" begin
            t1 = Trace(0, 1, [0, 1, 2])
            t1.meta.a = "A"
            t2 = Trace(0, 1, [1, 2, 3])
            t2.meta.a = "B"
            t2.meta.b = 2
            t3 = Trace(0, 1, [2, 3, 4])
            t3.meta.a = "C"
            t3.meta.b = 3
            t3.meta.c = 'c'
            t_merged = merge(t1, t2, t3)

            @testset "Out of place" begin
                @test t_merged.meta.a == "A"
                @test t_merged.meta.b == 2
                @test t_merged.meta.c == 'c'
                @test t1.meta.a == "A"
                @test t1.meta.b === missing
                @test t1.meta.c === missing
            end

            @testset "Empty traces" begin
                ts = [t1, t2, t3]
                @testset "Empty trace $i" for i in eachindex(ts)
                    ts′ = deepcopy(ts)
                    empty!(trace(ts′[i]))
                    @test merge(ts′).meta == Seis.SeisDict{Symbol,Any}(
                        :a=>"A", :b=>2, :c=>'c')
                end
            end

            @testset "In-place" begin
                merge!(t1, t2, t3)
                @test t1.meta == t_merged.meta
            end
        end

        @testset "Merge helpers" begin
            @testset "_find_gap_sections" begin
                @test Seis._find_gap_sections([]) == []
                @test Seis._find_gap_sections([0, 0, 0]) == [(1, 3)]
                @test Seis._find_gap_sections([0, 1, 2, 3, 0, 0]) == [(1, 1), (5, 6)]
                @test Seis._find_gap_sections([0, 1, 2, 3, 0, 0, 2]) == [(1, 1), (5, 6)]
                @test Seis._find_gap_sections([0, 0, 1, 2, 3, 0, 0, 2, 0]) == [(1, 2), (6, 7), (9, 9)]
            end

            @testset "_traces_are_sequential" begin
                @test Seis._traces_are_sequential(
                    [Trace(0, 1, [1, 2, 3]), Trace(3, 1, [4, 5])], [1, 2], 0.1)
                @test !Seis._traces_are_sequential(
                    [Trace(0, 1, [1, 2, 3]), Trace(3.5, 1, [4, 5])], [1, 2], 0.1)
                @test !Seis._traces_are_sequential(
                    [Trace(0, 1, [1, 2, 3]), Trace(4, 1, [4, 5])], [1, 2], 0.1)
                @test !Seis._traces_are_sequential(
                    [Trace(0, 1, [1, 2, 3]), Trace(2, 1, [4, 5])], [1, 2], 0.1)
                @test Seis._traces_are_sequential(
                    [Trace(3, 1, [4, 5]), Trace(0, 1, [1, 2, 3])], [2, 1], 0.1)
                @test_throws ArgumentError Seis._traces_are_sequential([Trace(0, 1, [])], [1], 0.1)
                @test_throws ArgumentError Seis._traces_are_sequential([], [], 0.1)
            end

            @testset "_channel_codes_are_equal" begin
                t1 = Trace(0, 1, [])
                t1.sta.net = "A"
                t1.sta.sta = "B"
                t1.sta.loc = "C"
                t1.sta.cha = "D"

                t2 = deepcopy(t1)

                @test Seis._channel_codes_are_equal(t1, t2)
                @test Seis._channel_codes_are_equal(t1, [t1, t2])

                t2.sta.loc = missing
                @test !Seis._channel_codes_are_equal(t1, t2)
                @test !Seis._channel_codes_are_equal(t1, [t1, t2])

            end
        end
    end
end
