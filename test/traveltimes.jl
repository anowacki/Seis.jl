using Test
using Seis

@testset "Picks" begin
    let t = Trace(0, 1, rand(2)), pick_time = rand(), pick_name = "pPKiKPPKiKP",
            pick_key = Symbol(pick_name)
        @test picks(t) == []
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
    end
end
