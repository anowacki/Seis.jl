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

    @testset "Position" begin
        let lon = rand(), lat = rand(), dep = rand(), (x, y, z) = rand(3)
            local g, c
            # Geographic
            g = Seis.Geographic(lon, lat, dep)
            @test g isa Seis.Geographic{Float64}
            @test typeof(g) <: Seis.Position
            @test typeof(g) <: Seis.Position{Float64}
            @test g == Seis.Geographic{Float64}(lon, lat, dep)
            @test g.lon == lon
            @test g.lat == lat
            @test g.dep == dep
            @test Seis.Geographic() == Seis.Geographic{Float64}(missing, missing, missing)
            @test Seis.Geographic(lon=lon, lat=lat, dep=dep) == Seis.Geographic{Float64}(lon, lat, dep)
            @test Seis.Geographic(lon=lon, lat=lat).dep === missing
            # Cartesian
            c = Seis.Cartesian(x, y, z)
            @test c isa Seis.Cartesian{Float64}
            @test typeof(c) <: Seis.Position
            @test typeof(c) <: Seis.Position{Float64}
            @test c == Seis.Cartesian{Float64}(x, y, z)
            @test c.x == x
            @test c.y == y
            @test c.z == z
            @test Seis.Cartesian() == Seis.Cartesian{Float64}(missing, missing, missing)
            @test Seis.Cartesian(x=x, y=y, z=z) == Seis.Cartesian{Float64}(x, y, z)
            @test Seis.Cartesian(x=x, y=y).z === missing
            # Arrays
            @test [g, g, g].lon == [lon, lon, lon]
            @test [c, c, c].x == [x, x, x]
            # Integer indexing
            @test g[1] == lon
            @test g[2] == lat
            @test g[3] == dep
            @test c[1] == x
            @test c[2] == y
            @test c[3] == z
        end
    end

    @testset "Event" begin
        # Geographic events
        let dt = now(), lon = rand(), lat = rand(), dep = rand(), id = "id",
                meta = Dict()
            local e, es
            @test (e = Event()) isa Event{Float64}
            @test all(ismissing, (e.lon, e.lat, e.dep, e.time))
            @test_throws MethodError Event("x")
            @test_throws TypeError Event{Int,Seis.Geographic{Int}}()
            e = Event(lon=lon, lat=lat, dep=dep, time=dt, id=id, meta=meta)
            @test e.pos == Seis.Geographic(lon, lat, dep)
            @test e.lon == lon
            @test e.lat == lat
            @test e.dep == dep
            @test e.time == dt
            @test e.id == id
            @test e.meta == meta
            e.meta.key_name = dt
            @test e.meta.key_name == dt
            # Arrays
            es = [deepcopy(e) for _ in 1:3]
            @test typeof(es) == Vector{Event{Float64, Seis.Geographic{Float64}}}
            @test length(es.lon) == 3
            @test es.lat == [lat, lat, lat]
            es.id = ["1", "2", "3"]
            @test es.id == ["1", "2", "3"]
            es.id = "A"
            @test es.id == ["A", "A", "A"]
            es .= e
            @test all(x->x==e, es)
        end

        # Cartesian events
        let dt = now(), x = rand(), y = rand(), z = rand(), id = "id",
                meta = Dict()
            local e, es
            @test CartEvent() isa Event{Float64, Seis.Cartesian{Float64}}
            @test CartEvent{Float32}() isa Event{Float32, Seis.Cartesian{Float32}}
            e = CartEvent()
            @test all(ismissing, (e.x, e.y, e.z, e.time, e.id))
            e = CartEvent(x=x, y=y, z=z, time=dt, id=id, meta=meta)
            @test e.pos == Seis.Cartesian(x, y, z)
            @test e.x == x
            @test e.y == y
            @test e.z == z
            @test e.time == dt
            @test e.id == id
            @test e.meta == meta
            e.meta.key_name = dt
            @test e.meta.key_name == dt
            # Arrays
            es = [deepcopy(e) for _ in 1:3]
            @test typeof(es) == Vector{Event{Float64, Seis.Cartesian{Float64}}}
            @test length(es.x) == 3
            @test es.y == [y, y, y]
            es.id = ["1", "2", "3"]
            @test es.id == ["1", "2", "3"]
            es.id = "A"
            @test es.id == ["A", "A", "A"]
            es .= e
            @test all(x->x===e, es)
        end
    end

    @testset "Station" begin
        # Geographic events
        let lon = rand(), lat = rand(), elev = rand(), net = "A", sta = "B",
                cha = "SX1", loc = "00", dep = rand(), azi = rand(), inc = rand(),
                meta = Dict(), dt = now()
            local s
            @test (s = Station()) isa Station{Float64,Seis.Geographic{Float64}}
            @test all(ismissing, (s.net, s.sta, s.loc, s.cha, s.lon, s.lat,
                                  s.dep, s.elev, s.azi, s.inc))
            @test_throws MethodError Station(1)
            s = Station(; net=net, sta=sta, loc=loc, cha=cha, lon=lon, lat=lat,
                dep=dep, elev=elev, azi=azi, inc=inc, meta=meta)
            @test s.net == net
            @test s.sta == sta
            @test s.loc == loc
            @test s.cha == cha
            @test s.lon == lon
            @test s.lat == lat
            @test s.dep == dep
            @test s.pos == Seis.Geographic(lon, lat, dep)
            @test s.elev == elev
            @test s.azi == azi
            @test s.inc == inc
            @test s.meta == meta
            s.meta.key_name = dt
            @test s.meta.key_name == dt
            s.azi *= 2
            @test s.azi == azi*2
            # Arrays
            ss = [deepcopy(s) for _ in 1:3]
            @test typeof(ss) == Vector{Station{Float64,Seis.Geographic{Float64}}}
            @test length(ss.lon) == 3
            @test ss.lat == [lat, lat, lat]
            ss.loc = ["1", "2", "3"]
            @test ss.loc == ["1", "2", "3"]
            ss.loc = "A"
            @test ss.loc == ["A", "A", "A"]
            ss .= s
            @test all(x->x===s, ss)
        end

        # Cartesian stations
        let x = rand(), y = rand(), z = rand(), elev = rand(), net = "A", sta = "B",
                cha = "SX1", loc = "00", dep = rand(), azi = rand(), inc = rand(),
                meta = Dict(), dt = now()
            local s
            @test (s = CartStation()) isa
                Station{Float64,Seis.Cartesian{Float64}}
            @test CartStation{Float32}() isa
                Station{Float32,Seis.Cartesian{Float32}}
            @test all(ismissing, (s.net, s.sta, s.loc, s.cha, s.x, s.y, s.z,
                                  s.elev, s.azi, s.inc))
            @test_throws MethodError CartStation(1)
            s = CartStation(; net=net, sta=sta, loc=loc, cha=cha, x=x, y=y, z=z,
                elev=elev, azi=azi, inc=inc, meta=meta)
            @test s.net == net
            @test s.sta == sta
            @test s.loc == loc
            @test s.cha == cha
            @test s.x == x
            @test s.y == y
            @test s.z == z
            @test s.pos == Seis.Cartesian(x, y, z)
            @test s.elev == elev
            @test s.azi == azi
            @test s.inc == inc
            @test s.meta == meta
            s.meta.key_name = dt
            @test s.meta.key_name == dt
            s.azi *= 2
            @test s.azi == azi*2
            # Arrays
            ss = [deepcopy(s) for _ in 1:3]
            @test typeof(ss) == Vector{Station{Float64,Seis.Cartesian{Float64}}}
            @test length(ss.x) == 3
            @test ss.x == [x, x, x]
            ss.loc = ["1", "2", "3"]
            @test ss.loc == ["1", "2", "3"]
            ss.loc = "A"
            @test ss.loc == ["A", "A", "A"]
            ss .= s
            @test all(x->x===s, ss)
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
            @test_throws Exception Trace(b, delta, -n)
        end

        let b = 1, delta = 1, tr = rand(1:3, 100), t = Trace(b, delta, tr)
            @test t isa Trace
            @test typeof(t) <: Seis.AbstractTrace
            @test typeof(t) <: Trace
            @test typeof(t) == Trace{Float64,Vector{Float64},Seis.Geographic{Float64}}
        end

        let t = Trace{Float32,Vector{Float32}}(rand(), rand(), rand(3))
            @test t isa Trace{Float32,Vector{Float32},Seis.Geographic{Float32}}
        end
        @test Trace{Float32}(0, 1, 100) isa Trace{Float32,Vector{Float32},Seis.Geographic{Float32}}

        @test_throws TypeError Trace{Int,Vector{Int}}(1, 1, rand(1:3, 3))

        # Broadcasting traces as scalars
        @test nsamples.(Trace(0, 1, rand(2))) == 2

        # Arrays
        let b = rand(), t = Trace(b, rand(), rand(3))
            ts = [deepcopy(t) for _ in 1:3]
            @test ts isa Array{<:Trace,1}
            @test length(ts.b) == 3
            @test ts.b == [b, b, b]
            ts.delta = [1, 2, 3]
            @test ts.delta == [1, 2, 3]
            ts.delta = 2
            @test ts.delta == [2, 2, 2]
            @test ts.evt isa Vector{Event{Float64,Seis.Geographic{Float64}}}
            ts′ = deepcopy(ts)
            for tt in ts
                tt.evt.id = "A"
            end
            @test ts.evt.id == ["A", "A", "A"]
            ts′.evt.id = "A"
            @test ts == ts′
            ts.evt = t.evt
            @test all(x->x===t.evt, ts.evt)
            ts .= t
            @test all(x->x===t, ts)
        end
    end
end
