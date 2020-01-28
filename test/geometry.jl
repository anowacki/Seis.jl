using Test
using Dates
using Seis

@testset "Geometry" begin
    # Geographic
    e = Event(lon=0, lat=0)
    s = Station(lon=0, lat=45)
    t = Trace(0, 1, 1)
    t.sta = s
    t.evt = e
    @test azimuth(e, s) ≈ 0 atol=√eps(Float64)
    @test azimuth(e, s) == azimuth(t)
    @test backazimuth(s, e) ≈ 180
    @test backazimuth(s, e) == backazimuth(t)
    @test distance_deg(e, s) == distance_deg(t)
    @test distance_deg(e, s, sphere=true) == 45.0
    @test distance_km(e, s) == distance_km(t)
    @test distance_km(e, s, sphere=true) ≈ π/4*6371 rtol=0.01
    e.lon = missing
    @test_throws ArgumentError azimuth(e, s)
    @test_throws ArgumentError backazimuth(s, e)
    @test_throws MethodError incidence(e, s)
    @test_throws ArgumentError azimuth(t)
    @test_throws ArgumentError backazimuth(t)
    @test_throws MethodError incidence(t)
    # Cartesian
    t = CartTrace(0, 1, 1)
    s = t.sta
    e = t.evt
    s.x, s.y, s.z = -1, 1, 0
    e.x, e.y, e.z = 0, 0, -1
    @test azimuth(e, s) ≈ 315
    @test azimuth(e, s) == azimuth(t)
    @test backazimuth(s, e) ≈ 135
    @test backazimuth(s, e) == backazimuth(t)
    @test_throws MethodError distance_deg(e, s)
    @test_throws MethodError distance_deg(t)
    @test distance_km(e, s) ≈ √2/1000
    @test distance_km(e, s) == distance_km(t)
    @test incidence(e, s) ≈ rad2deg(atan(√2))
    e.x = missing
    @test_throws ArgumentError azimuth(e, s)
    @test_throws ArgumentError backazimuth(s, e)
    @test_throws ArgumentError incidence(e, s)
    @test_throws ArgumentError azimuth(t)
    @test_throws ArgumentError backazimuth(t)
    @test_throws ArgumentError incidence(t)
end
