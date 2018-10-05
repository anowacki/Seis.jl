using Test
using Seis
import TauPy

@testset "Picks" begin
    let t = Trace(0, 1, rand(2)), pick_time = rand(), pick_name = "pPKiKPPKiKP"
        @test picks(t) == []
        add_pick!(t, pick_time, pick_name)
        @test length(picks(t)) == 1
        p = picks(t)
        @test p[1][1] == pick_time
        @test p[1][2] == pick_name
        add_pick!(t, pick_time+1)
        @test length(picks(t)) == 2
        p = picks(t)
        @test p[2][1] == pick_time+1
        @test ismissing(p[2][2])
        @test_throws ArgumentError add_picks!(t, "S")
        t.evt.lon, t.evt.lat, t.evt.dep = 0, 0, 0
        t.sta.lon, t.sta.lat = 10, 10
        # Adding travel time picks from TauPy
        add_picks!(t, "Sn", model="ak135")
        @test picks(t)[end][1] ≈ 358.4306 atol=0.0001
        @test picks(t)[end][2] == "Sn"
    end
end

@testset "Arrivals" begin
    let p = TauPy.travel_time(100, 100)

    end
end
