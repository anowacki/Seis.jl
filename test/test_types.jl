using Test
using Dates
using Seis

@testset "Types" begin
    @testset "SeisDict" begin
        let
            local d, dd
            d = Seis.SeisDict(:a=>1)
            @test d[:a] == 1
            @test d.a == 1
            @test d.b === missing
            d.b = 2
            @test d[:b] == d.b == 2
            d.a = d.b = missing
            @test haskey(d, :a) == false
            @test haskey(d, :b) == false
            @test collect(keys(d)) == []
            dd = [Seis.SeisDict(:a=>i) for i in 1:3]
            @test dd.a == [1,2,3]
            dd.a = [2,3,4]
            @test dd.a == [2,3,4]
            @test dd[1].a == 2
            dd.a = "x"
            @test dd.a == ["x", "x", "x"]
            dd.a = missing
            @test collect.(keys.(dd)) == [[], [], []]
        end
    end

    @testset "Event" begin
        let dt = now(), lon = rand(), lat = rand(), dep = rand(), id = "id",
                dt = now(), meta = Dict()
            local e
            @test (e = Event()) isa Event{Float64,String}
            @test all(ismissing, (e.lon, e.lat, e.dep, e.time))
            @test_throws MethodError Event("x")
            @test_throws TypeError Event{Int,String}()
            e = Event(lon, lat, dep, dt, id, meta)
            @test e.lon == lon
            @test e.lat == lat
            @test e.dep == dep
            @test e.time == dt
            @test e.id == id
            @test e.meta == meta
            e.meta.key_name = dt
            @test e.meta.key_name == dt
        end
    end

    @testset "Station" begin
        let lon = rand(), lat = rand(), elev = rand(), net = "A", sta = "B",
                cha = "SX1", loc = "00", dep = rand(), azi = rand(), inc = rand(),
                meta = Dict(), dt = now()
            local s
            @test (s = Station()) isa Station{Float64,String}
            @test all(ismissing, (s.net, s.sta, s.loc, s.cha, s.dep, s.elev, s.azi, s.inc))
            @test_throws MethodError Station(1)
            s = Station(net, sta, loc, cha, lon, lat, dep, elev, azi, inc, meta)
            @test s.net == net
            @test s.sta == sta
            @test s.loc == loc
            @test s.cha == cha
            @test s.lon == lon
            @test s.lat == lat
            @test s.dep == dep
            @test s.elev == elev
            @test s.azi == azi
            @test s.inc == inc
            @test s.meta == meta
            s.meta.key_name = dt
            @test s.meta.key_name == dt
            s.azi *= 2
            @test s.azi == azi*2
        end
    end

    @testset "Trace" begin
        let b = rand(), delta = rand(), n = rand(1:100), tr = rand(n),
                t = Trace(b, delta, tr)
            @test nsamples(t) == n
            @test all(tr .== trace(t))
            @test t.b == b
            @test t.delta == delta

            @test_throws MethodError Trace()
            @test_throws MethodError Trace(b)
            @test_throws MethodError Trace(b, delta)
            @test_throws ArgumentError Trace(b, -delta, tr)
            @test_throws MethodError Trace(b, delta, b)
        end

        let b = rand(), delta = rand(), n = rand(1:100), t = Trace(b, delta, n)
            @test t.b == b
            @test t.delta == delta
            @test nsamples(t) == n
            @test_throws ErrorException Trace(b, delta, -n)
        end

        let b = 1, delta = 1, tr = rand(1:3, 100), t = Trace(b, delta, tr)
            @test t isa Trace
            @test typeof(t) <: Seis.AbstractTrace
            @test typeof(t) <: Trace
            @test typeof(t) == Trace{Float64,Vector{Float64},String}
        end

        let t = Trace{Float32,Vector{Float32},String}(rand(), rand(), rand(3))
            @test t isa Trace{Float32,Vector{Float32},String}
        end
        @test Trace{Float32,Vector{Float16},SubString}(0, 1, 10) isa Trace{Float32,Vector{Float16},SubString}
        @test Trace{Float32}(0, 1, 100) isa Trace{Float32,Vector{Float32},String}

        @test_throws TypeError Trace{Int,Vector{Int},String}(1, 1, rand(1:3, 3))
    end
end
