using Test
using Seis

@testset "Picks" begin
    @testset "Construction" begin
        @testset "$T" for T in (Float32, Float64)
            @test Seis.Pick{T}(2) isa Seis.Pick{T}
            @test Seis.Pick{T}(1).time === T(1)
            @test Seis.Pick{T}(1).name === missing
            @test Seis.Pick{T}(1, "A").time === T(1)
            @test Seis.Pick{T}(1, "A").name === "A"
            @test Seis.Pick{T}(; time=1) == Seis.Pick{T}(1, missing)
            @test Seis.Pick{T}(; time=2, name="B") == Seis.Pick{T}(2, "B")
        end
        @testset "Untyped" begin
            default_float = Float64
            @test Seis.Pick(1) isa Seis.Pick{default_float}
            @test Seis.Pick(1, "A") isa Seis.Pick{default_float}
            @test Seis.Pick(1).time === default_float(1)
            @test Seis.Pick(1).name === missing
            @test Seis.Pick(1, "A").time === default_float(1)
            @test Seis.Pick(1, "A").name ==="A"
            @test Seis.Pick(; time=-1) == Seis.Pick(-1)
            @test Seis.Pick(; name="A", time=2) == Seis.Pick(2, "A")
        end
    end

    @testset "Conversion" begin
        @testset "$T" for T in (Float32, Float64)
            let time = rand(T), name = "XYZ"
                @test convert(Seis.Pick{T}, time) == Seis.Pick{T}(time)
                @test convert(Seis.Pick{T}, (time, name)) == Seis.Pick{T}(time, name)
                @test convert(Seis.Pick{T}, (time=time, name=name)) == Seis.Pick{T}(time, name)
                from_type = T == Float32 ? Float64 : Float32
                @test convert(Seis.Pick{T}, Seis.Pick{from_type}(1, "A")) isa Seis.Pick{T}
                @test convert(Seis.Pick{T}, Seis.Pick{from_type}(1, "A")).time === T(1)
            end
        end
        @testset "Untyped" begin
            @test convert(Seis.Pick, 1) == Seis.Pick{Float64}(1)
            @test convert(Seis.Pick, (1, "A")) == Seis.Pick{Float64}(1, "A")
            @test convert(Seis.Pick, Seis.Pick{Float32}(1, "A")) == Seis.Pick{Float64}(1, "A")
            @test convert(Seis.Pick, Seis.Pick(1, "A")) === Seis.Pick(1, "A")
        end
    end

    # Picks can be iterated to give their time, then name
    @testset "Iteration" begin
        let time = rand(), name = "XYZ"
            p = Seis.Pick(time, name)
            v, i = iterate(p)
            @test v == time
            v, i = iterate(p, i)
            @test v == name
            @test iterate(p, i) === nothing
            @test first(p) == time
            @test last(p) == name
            a, b = p
            @test a == time
            @test b == name
        end
        let time = rand()
            p = Seis.Pick(time)
            a, b = p
            @test (a, b) === (time, missing)
            p′ = Seis.Pick(time, "A")
            a, b = p′
            @test (a, b) === (time, "A")
        end
    end

    # Picks are broadcast as single objects, not an iterable collection
    @testset "Broadcasting" begin
        let a = Vector{Any}(undef, 3), p = [Seis.Pick(rand()) for _ in 1:3]
            a .= p[1]
            @test all(a .== p[1])
            a .= p
            @test all(a .== p)
        end
    end

    @testset "Comparison" begin
        @test Seis.Pick(2) == Seis.Pick(2) === Seis.Pick(2)
        @test Seis.Pick(2, "A") != Seis.Pick(2)
        @test Seis.Pick(3, "A") != Seis.Pick(3, "B")
        @test Seis.Pick(4, "C") === Seis.Pick(4, "C")
    end

    let t = Trace(0, 1, rand(2)), pick_time = rand(), pick_name = "pPKiKPPKiKP",
            pick_key = Symbol(pick_name)
        add_pick!(t, pick_time, pick_name)
        @test length(picks(t)) == 1
        @test t.picks[pick_key].time == pick_time
        @test t.picks[pick_key].name == pick_name
        @test t.picks[pick_key] == t.picks.pPKiKPPKiKP
        add_pick!(t, pick_time+1)
        @test t.picks[1].time == pick_time + 1
        @test ismissing(t.picks[1].name)

        # Missing is returned for key not present
        @test !haskey(t.picks, 2)
        @test !haskey(t.picks, :A)
        @test ismissing(t.picks[2])
        @test ismissing(t.picks.A)

        # Setting an entry to missing removes it
        t.picks.A = 1, "B"
        t.picks.B = 2, "C"
        t.picks[2] = 3
        @test t.picks[:A] == t.picks.A == Seis.Pick(1, "B")
        @test t.picks[2] == Seis.Pick(3)
        t.picks.A = missing
        t.picks[:B] = missing
        t.picks[2] = missing
        @test !haskey(t.picks, :A)
        @test ismissing(t.picks.A)
        @test !haskey(t.picks, :B)
        @test ismissing(t.picks.B)
        @test !haskey(t.picks, 2)
        @test ismissing(t.picks[2])

        # Removing picks
        clear_picks!(t)
        @test length(picks(t)) == 0

        # Picks ordered by time by default
        add_pick!(t, pick_time, pick_name)
        p = picks(t)
        @test all(p .=== picks(t, sort=:time)) # all(.===) ∵ have missings
        @test p[1][1] == p[1].time == pick_time
        @test p[1][2] == p[1].name == pick_name
        add_pick!(t, pick_time+1)
        @test length(picks(t)) == 2
        p = picks(t)
        @test p[2].time == pick_time+1
        @test ismissing(p[2].name)

        # Picks sorted by name
        clear_picks!(t)
        add_pick!.(t, (1, 2), ("B", "A"))
        add_pick!(t, 3)
        @test all(picks(t, sort=:name).name .=== [missing, "A", "B"])

        # Removing all picks
        clear_picks!(t)
        @test length(t.picks) == length(picks(t)) == 0
        @test t.picks[pick_key] === missing

        # Correct return type for empty sets of picks
        @testset "picks no key" begin
            @testset "Eltype $T" for T in (Float32, Float64)
                let t = Trace{T}(0, 1, 10)
                    @test isempty(picks(t, :nokey))
                    @test picks(t, :nokey) isa Vector{Seis.Pick{T}}
                    @test isempty(picks(t, "noname"))
                    @test picks(t, "noname") isa Vector{Seis.Pick{T}}
                    @test isempty(picks(t, r"nomatch"))
                    @test picks(t, r"nomatch") isa Vector{Seis.Pick{T}}
                end
            end
        end
    end

    @testset "add_pick!" begin
        @testset "Relative time" begin
            let t = Trace(0, 1, 1), p = Seis.Pick(1, "A")
                add_pick!(t, p)
                @test haskey(t.picks, :A)
                @test t.picks.A == p
                add_pick!(t, p)
                @test haskey(t.picks, :A_1)
                @test t.picks.A_1 == p

                add_pick!(t, p, "B")
                @test haskey(t.picks, :B)
                @test t.picks.B == Seis.Pick(p.time, "B")
                add_pick!(t, p, missing)
                @test haskey(t.picks, 1)
                @test t.picks[1] == Seis.Pick(p.time)

                add_pick!(t, p, missing)
                @test haskey(t.picks, 2)
                @test t.picks[2] == t.picks[1]
            end
        end

        @testset "DateTime" begin
            t = Trace(0, 1, 0)
            t.evt.time = DateTime(3000)
            @test add_pick!(t, DateTime(3000, 1, 1, 0, 0, 1, 360)) == Seis.Pick(time=1.36, name=missing)
            @test haskey(t.picks, 1)
            @test t.picks[1] == Seis.Pick(time=1.36, name=missing)

            @test add_pick!(t, DateTime(2999, 12, 31, 23, 59, 50), "B") == Seis.Pick(time=-10, name="B")
            @test t.picks[:B] == Seis.Pick(time=-10.0, name="B")

            @test picks(t) == [Seis.Pick(-10.0, "B"), Seis.Pick(1.36)]
        end
    end
end
