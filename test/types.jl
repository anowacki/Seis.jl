using Test
using Dates
using Seis

@static if VERSION < v"1.2"
    hasproperty(d, k) = k in propertynames(d)
end

@testset "Types" begin
    @testset "SeisDict" begin
        @testset "Single dicts" begin
            d = Seis.SeisDict(:a=>1)
            @test d[:a] == 1
            @test haskey(d, :a) == true
            @test hasproperty(d, :a) == true
            @test propertynames(d) == [:a]
            @test d.a == 1
            @test d.b === missing
            @test hasproperty(d, :b) == false
            d.b = 2
            @test d[:b] == d.b == 2
            @test propertynames(d) == [:a, :b]
            d.a = d.b = missing
            @test haskey(d, :a) == false
            @test haskey(d, :b) == false
            @test collect(keys(d)) == []
            @test hasproperty(d, :a) == false
        end
        
        @testset "Arrays" begin
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
        @testset "Geographic" begin
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
            @testset "meta arrays" begin
                es = [Event() for _ in 1:3]
                @test es.meta == [Seis.SeisDict{Symbol,Any}() for _ in 1:3]
                es.meta.x = [1, 2, 3]
                @test es[1].meta.x == 1
                @test es.meta.x == [1, 2, 3]
                @test all(x -> x === missing, es.meta.y)
                es.meta.x = missing
                @test all(x -> x === missing, es.meta.x)
            end
        end

        @testset "Cartesian" begin
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
            @testset "meta arrays" begin
                es = [CartEvent() for _ in 1:3]
                @test es.meta == [Seis.SeisDict{Symbol,Any}() for _ in 1:3]
                es.meta.x = [1, 2, 3]
                @test es[1].meta.x == 1
                @test es.meta.x == [1, 2, 3]
                @test all(x -> x === missing, es.meta.y)
                es.meta.x = missing
                @test all(x -> x === missing, es.meta.x)
            end
        end
    end

    @testset "Station" begin
        @testset "Geographic" begin
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
            @testset "meta arrays" begin
                ss = [Station() for _ in 1:3]
                @test ss.meta == [Seis.SeisDict{Symbol,Any}() for _ in 1:3]
                ss.meta.x = [1, 2, 3]
                @test ss[1].meta.x == 1
                @test ss.meta.x == [1, 2, 3]
                @test all(x -> x === missing, ss.meta.y)
                ss.meta.x = missing
                @test all(x -> x === missing, ss.meta.x)
            end
        end

        @testset "Cartesian" begin
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
            @testset "meta arrays" begin
                ss = [CartStation() for _ in 1:3]
                @test ss.meta == [Seis.SeisDict{Symbol,Any}() for _ in 1:3]
                ss.meta.x = [1, 2, 3]
                @test ss[1].meta.x == 1
                @test ss.meta.x == [1, 2, 3]
                @test all(x -> x === missing, ss.meta.y)
                ss.meta.x = missing
                @test all(x -> x === missing, ss.meta.x)
            end
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
            @test all(x->x==t, ts)
            # Test for aliasing
            trace(t) .= -1
            @test all(x->all(trace(x) .== -1), ts)
        end
    end

    # Conversion tests for Picks are in traveltimes.jl
    # TODO: Move Pick tests here.
    @testset "Conversion" begin
        @testset "$Geom" for Geom in (Seis.Geographic, Seis.Cartesian)
            Ts = (Float32, Float64)
            @testset "$Tin -> $Tout" for Tin in Ts, Tout in Ts
                @testset "Position $Geom" begin
                    intype = Geom{Tin}
                    outtype = Geom{Tout}
                    @test convert(outtype, intype(1, 2, 3)) isa outtype
                    @test getfield(convert(outtype, intype(1, 2, 3)), 1) == 1
                    @test convert(outtype, intype(missing, missing, missing)) isa outtype
                    @test getfield(convert(outtype, intype(missing, 2, 3)), 1) === missing
                    # Identity 'conversions' give the same object
                    if intype == outtype
                        in = intype(1, 2, 3)
                        @test convert(outtype, in) === in
                    end
                end

                # Events and Stations
                @testset "$ES" for ES in (Event, Station)
                    intype = ES{Tin,Geom{Tin}}
                    outtype = ES{Tout,Geom{Tout}}
                    kwargs = Geom == Seis.Geographic ? (lon=1, lat=2, dep=3) : (x=1, y=2, z=3)
                    if ES == Station
                        kwargs = (net="A", sta="B", loc="C", cha="D", elev=-1,
                                  azi=rand(), inc=rand(), meta=Dict(:x=>"X"), kwargs...)
                    else
                        kwargs = (id="XYZ", time=now(), meta=Dict(:y=>"Y"), kwargs...)
                    end
                    @test convert(outtype, intype()) isa outtype
                    @test convert(outtype, intype(; kwargs...)) isa outtype
                    posfield = ES == Event ? 1 : 5
                    @test getfield(convert(outtype, intype(; kwargs...)), posfield) == Geom{Tout}(1, 2, 3)
                    # Identity 'conversions' give the same object
                    if intype == outtype
                        in = intype(; kwargs...)
                        @test convert(outtype, in) === in
                    end
                end

                @testset "Trace $Geom" begin
                    Vs = [Vector{TT} for TT in Ts]
                    @testset "Data type $Vin -> $Vout" for Vin in Vs, Vout in Vs
                        intype = Trace{Tin, Vin, Geom{Tin}}
                        outtype = Trace{Tout, Vout, Geom{Tout}}
                        vtype = eltype(Vout)
                        data = vtype[1, 2, 3]
                        tin = intype(0, 1, data)
                        tout = convert(outtype, tin)
                        @test tout isa outtype
                        @test trace(tout) isa Vout
                        @test all(trace(tin) .≈ trace(tout))
                        # Check aliasing is as expected
                        trace(tin) .= -1
                        if Vin == Vout
                            @test trace(tin) === trace(tout)
                            if intype == outtype
                                @test tin === tout
                            else
                                @test tin !== tout
                            end
                        else
                            @test trace(tin) !== trace(tout)
                        end
                    end
                end
            end
        end
    end

    # Test that views of arrays of our types do not cause errors
    @testset "Views" begin
        @testset "$T" for T in (Event, Station)
            es = [T(lon=rand()) for _ in 1:3]
            ev = view(es, 2:3)
            @test ev.lon == es[2:3].lon
            @test ev[1] == es[2]
        end

        @testset "Pick" begin
            ps = [Seis.Pick(rand(), name) for name in ("A", "B", "C")]
            pv = view(ps, 2:3)
            @test pv.time == ps[2:3].time
            @test pv[1] == ps[2]
        end

        @testset "Trace" begin
            ts = [Trace(rand(), 1, 0) for _ in 1:3]
            tv = view(ts, 2:3)
            @test tv.b == ts[2:3].b
            @test tv[1] == ts[2]
        end
    end

    # Hashing, isequal and uniqueness
    function test_hash_isequal(x, x′, y)
        @test isequal(x, x′)
        @test !isequal(x, y)
        @test hash(x) == hash(x′)
        @test hash(x) != hash(y)
        @test unique([x, x′, y]) == [x, y]
    end

    @testset "Hash and isequal" begin
        @testset "Cartesian" begin
            c = [Seis.Cartesian(x=1) for _ in 1:3]
            c[3].y = 2
            test_hash_isequal(c...)
        end
        @testset "Geographic" begin
            g = [Seis.Geographic(lon=1) for _ in 1:3]
            g[3].lat = 2
            test_hash_isequal(g...)
        end
        @testset "$T" for T in (Event, Station)
            e = [T(lon=1) for _ in 1:3]
            e[3].meta.a = 1
            test_hash_isequal(e...)
            e[3].meta.a = missing
            e[3].dep = 2
            test_hash_isequal(e...)
        end
        @testset "Trace" begin
            t = [Trace(0, 1, 0) for _ in 1:3]
            push!(trace(t[3]), 1)
            test_hash_isequal(t...)
            empty!(trace(t[3]))
            t[3].b = 1
            test_hash_isequal(t...)
            t[3].b = t[1].b
            t[3].sta.lon = 1
            test_hash_isequal(t...)
            t[3].sta.lon = missing
            t[3].evt.time = now()
            test_hash_isequal(t...)
            t[3].evt.time = missing
            t[3].meta.a = 1
            test_hash_isequal(t...)
        end
    end
end
