# Test for interoperation with SeisIO

using Test
using Seis
import Seis.SeisIOIO
import SeisIO
using Dates: DateTime

@testset "SeisIOIO" begin
    @testset "Parsing" begin
        let x = rand(Float32, 10), chan = SeisIO.SeisChannel(
                    id="A.B.C.D", fs=10.0, t=[1 0; 10 0], x=x,
                    loc=SeisIO.GeoLoc(lon=1.0, lat=2.0, el=3.0, dep=4.0, az=5.0,
                                      inc=6.0, datum="WGS84"))
            # Default trace type
            @test SeisIOIO.seisio_date(chan) == DateTime(1970)
            @test length(SeisIOIO.parse_seisio(chan)) == 1
            t = first(SeisIOIO.parse_seisio(chan))
            @test t isa Seis.Trace{Float32, Vector{Float32}, Seis.Geographic{Float32}}
            @test t.sta.net == "A"
            @test t.sta.sta == "B"
            @test t.sta.loc == "C"
            @test t.sta.cha == "D"
            @test t.sta.lon == 1
            @test t.sta.lat == 2
            @test t.sta.elev == 3
            @test t.sta.dep == Float32(0.004)
            @test t.sta.azi == 5
            @test t.sta.inc == 6
            @test t.sta.meta.datum === missing
            @test nsamples(t) == 10
            @test t.delta == Float32(0.1)
            @test trace(t) == x
            chan.loc.datum = "XYZ"
            @test SeisIOIO.parse_seisio(chan)[1].sta.meta.datum == "XYZ"
            # Conversion to different geometry doesn't preserver location
            t′ = SeisIOIO.parse_seisio(CartTrace{Float64, Vector{Float64}}, chan)
            @test length(t′) == 1
            @test t′[1] isa CartTrace{Float64, Vector{Float64}}
            @test all(x -> x===missing, (t′[1].sta.x, t′[1].sta.y, t′[1].sta.z))
            # Cartesian conversion
            chan.loc = SeisIO.XYLoc(x=1.0, y=2.0, z=4.0, az=5.0, inc=6.0, ox=7.0,
                oy=8.0, oz=9.0, datum="XYZ")
            t″ = SeisIOIO.parse_seisio(CartTrace{Float32, Vector{Float32}}, chan)
            @test t″ isa Vector{CartTrace{Float32, Vector{Float32}}}
            sta = t″[1].sta
            @test sta.x == 1
            @test sta.y == 2
            @test sta.z == 4
            @test sta.azi == 5
            @test sta.inc == 6
            @test sta.meta.origin == (x=7.0, y=8.0, z=9.0)
            @test sta.meta.datum == "XYZ"
            # Missing fields
            chan.loc = SeisIO.GeoLoc()
            let sta = SeisIOIO.parse_seisio(chan)[1].sta
                @test all(x -> x===missing, (sta.lon, sta.lat, sta.elev, sta.dep, sta.inc, sta.azi))
            end
            chan.loc = SeisIO.XYLoc()
            let sta = SeisIOIO.parse_seisio(CartTrace{Float64, Vector{Float64}}, chan)[1].sta
                @test all(x -> x===missing, (sta.x, sta.y, sta.z, sta.azi, sta.inc))
            end
        end
    end

    @testset "Gaps" begin
        let chan = SeisIO.SeisChannel()
            chan.id = "AN.ABC..XYZ"
            chan.fs = 100
            chan.t = [1 1000; 101 1000; 151 -2000; 200 0]
            chan.x = rand(200)
            @test collect(SeisIOIO.gaps(chan)) == SeisIOIO.Gap.([101, 151], [1000, -2000])
            @test SeisIOIO.total_offset(chan) ≈ -0.001 atol = 1e-6
            @test SeisIOIO.largest_offset(chan) ≈ 0.002 atol = 1e-6
            # Allow all gaps
            t = SeisIOIO.parse_seisio(chan, maximum_gap=Inf, maximum_offset=Inf)
            @test length(t) == 1
            @test nsamples(t[1]) == 200
            @test t[1].meta.file === missing
            # Skip single-sample gaps that add up to a sample offset or less
            @test SeisIOIO.parse_seisio(chan) == t
            # Get all gaps
            t′ = SeisIOIO.parse_seisio(chan, maximum_gap=0)
            @test length(t′) == 3
            @test nsamples.(t′) == [100, 50, 50]
            @test sum(nsamples, t′) == sum(nsamples, t)
        end
    end
    @testset "IDs" begin
        let chan = SeisIO.SeisChannel(fs=20.0, x=rand(2), t=[1 0; 2 0])
            @test_logs (:warn, "channel code is not in an expected format") SeisIOIO.parse_seisio(chan)
            chan.id = "A.B.C.D"
            @test channel_code(SeisIOIO.parse_seisio(chan)[1]) == "A.B.C.D"
        end
    end
end
