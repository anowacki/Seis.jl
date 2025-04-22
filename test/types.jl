using Test
using Dates
using Seis
import FFTW

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
        let lon = rand(), lat = rand(), elev = rand(), (x, y, z) = rand(3)
            local g, c
            @testset "Geographic" begin
                g = Seis.Geographic(lon, lat, elev)
                @test g isa Seis.Geographic{Float64}
                @test typeof(g) <: Seis.Position
                @test typeof(g) <: Seis.Position{Float64}
                @test g == Seis.Geographic{Float64}(lon, lat, elev)
                @test g.lon == lon
                @test g.lat == lat
                @test g.elev == elev
                @test Seis.Geographic() == Seis.Geographic{Float64}(missing, missing, missing)
                @test Seis.Geographic(lon=lon, lat=lat, elev=elev) == Seis.Geographic{Float64}(lon, lat, elev)
                @test Seis.Geographic(lon=lon, lat=lat).elev === missing
                @test Seis.Geographic(dep=-0.001, elev=1) == Seis.Geographic(missing, missing, 1)
                @test_throws ArgumentError Seis.Geographic(dep=1, elev=1)

                @testset "setproperty" begin
                    g′ = deepcopy(g)
                    g′.lon = -lon
                    @test g′.lon == -lon
                end

                @testset "Depth conversion" begin
                    @test g.dep == -g.elev/1000
                    g′ = deepcopy(g)
                    g′.dep = 1
                    @test g′.elev == -1000
                end

                @testset "propertynames" begin
                    @test propertynames(Seis.Geographic()) ==
                        propertynames(Seis.Geographic(), false)
                    @test propertynames(Seis.Geographic()) == (:lon, :lat, :elev, :dep)
                end
            end

            @testset "Cartesian" begin
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

                @testset "setproperty" begin
                    c′ = deepcopy(c)
                    c′.y = -y
                    @test c′.y == -y
                end

                @testset "propertynames" begin
                    @test propertynames(Seis.Cartesian()) ==
                        propertynames(Seis.Cartesian(), false)
                    @test propertynames(Seis.Cartesian()) == (:x, :y, :z)
                end
            end

            @testset "Arrays" begin              
                @test [g, g, g].lon == [lon, lon, lon]
                @test [g, g].elev == [elev, elev]
                @test [g].dep == [-elev/1000]
                @test [c, c, c].x == [x, x, x]
            end

            @testset "Integer indexing" begin                
                @test g[1] == lon
                @test g[2] == lat
                @test g[3] == elev
                @test c[1] == x
                @test c[2] == y
                @test c[3] == z
            end
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
                @test e.pos == Seis.Geographic(; lon=lon, lat=lat, dep=dep)
                @test e.lon == lon
                @test e.lat == lat
                @test e.dep ≈ dep rtol=1e-10
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

            @testset "propertynames" begin
                @test propertynames(Event()) ==
                    propertynames(Event(), false)
                @test propertynames(Event()) == (:lon, :lat, :dep, :time, :id, :meta)
                @test propertynames(Event(), true) == (:lon, :lat, :dep, :pos, :time, :id, :meta)
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

            @testset "propertynames" begin
                @test propertynames(CartEvent()) ==
                    propertynames(CartEvent(), false)
                @test propertynames(CartEvent()) == (:x, :y, :z, :time, :id, :meta)
                @test propertynames(CartEvent(), true) == (:x, :y, :z, :pos, :time, :id, :meta)
            end
        end
    end

    @testset "Station" begin
        @testset "Geographic" begin
            let lon = rand(), lat = rand(), elev = rand(), net = "A", sta = "B",
                    cha = "SX1", loc = "00", azi = rand(), inc = rand(),
                    meta = Dict(), dt = now()
                local s
                @test (s = Station()) isa Station{Float64,Seis.Geographic{Float64}}
                @test all(ismissing, (s.net, s.sta, s.loc, s.cha, s.lon, s.lat,
                                      s.elev, s.azi, s.inc))
                @test_throws MethodError Station(1)
                s = Station(; net=net, sta=sta, loc=loc, cha=cha, lon=lon, lat=lat,
                    elev=elev, azi=azi, inc=inc, meta=meta)
                @test s.net == net
                @test s.sta == sta
                @test s.loc == loc
                @test s.cha == cha
                @test s.lon == lon
                @test s.lat == lat
                @test s.dep == -elev/1000
                @test s.pos == Seis.Geographic(lon, lat, s.elev)
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

            @testset "propertynames" begin
                @test propertynames(Station()) ==
                    propertynames(Station(), false)
                @test propertynames(Station()) ==
                    (:lon, :lat, :elev, :net, :sta, :loc, :cha, :azi, :inc, :meta)
                @test propertynames(Station(), true) ==
                    (:lon, :lat, :elev, :net, :sta, :loc, :cha, :pos, :azi, :inc, :meta)
            end
        end

        @testset "Cartesian" begin
            let x = rand(), y = rand(), z = rand(), net = "A", sta = "B",
                    cha = "SX1", loc = "00", azi = rand(), inc = rand(),
                    meta = Dict(), dt = now()
                local s
                @test (s = CartStation()) isa
                    Station{Float64,Seis.Cartesian{Float64}}
                @test CartStation{Float32}() isa
                    Station{Float32,Seis.Cartesian{Float32}}
                @test all(ismissing, (s.net, s.sta, s.loc, s.cha, s.x, s.y, s.z,
                                      s.azi, s.inc))
                @test_throws MethodError CartStation(1)
                s = CartStation(; net=net, sta=sta, loc=loc, cha=cha, x=x, y=y, z=z,
                    azi=azi, inc=inc, meta=meta)
                @test s.net == net
                @test s.sta == sta
                @test s.loc == loc
                @test s.cha == cha
                @test s.x == x
                @test s.y == y
                @test s.z == z
                @test s.pos == Seis.Cartesian(x, y, z)
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

            @testset "propertynames" begin
                @test propertynames(CartStation()) ==
                    propertynames(CartStation(), false)
                @test propertynames(CartStation()) ==
                    (:x, :y, :z, :net, :sta, :loc, :cha, :azi, :inc, :meta)
                @test propertynames(CartStation(), true) ==
                    (:x, :y, :z, :net, :sta, :loc, :cha, :pos, :azi, :inc, :meta)
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

            @test_throws MethodError Trace(b)
            @test_throws MethodError Trace(b, delta)
            @test_throws ArgumentError Trace(b, -delta, tr)
            @test_throws MethodError Trace(b, delta, b)

            @testset "Keyword constructor" begin
                @testset "$T $S $E" for (T, S, E) in (
                    (Trace, Station, Event),
                    (CartTrace, CartStation, CartEvent)
                )
                    @testset "Defaults" begin
                        @test T() == T(0, 1, [])
                    end

                    @testset "Specified kwargs" begin
                        sta = S(sta="ABC")
                        evt = E(id="123")
                        picks = Seis.SeisDict{Union{Int,Symbol},Seis.Pick{Seis.DEFAULT_FLOAT}}(
                            :a => Seis.Pick(1, "a"),
                            2 => Seis.Pick(2, "b")
                        )
                        meta = Seis.SeisDict(Dict{Symbol,Any}(:test=>0.1))
                        @test T(b=b).b == b
                        @test trace(T(data=tr)) == tr
                        @test nsamples(T(n=n-1)) == n - 1
                        @test T(delta=delta).delta == delta
                        @test T(evt=evt).evt == evt
                        @test T(sta=sta).sta == sta
                        @test T(meta=meta).meta == meta
                        @test T(picks=picks).picks == picks

                        @test (
                            @test_logs (
                                :warn, "ignoring keyword argument `n` as `data` was passed in"
                            ) T(n=10, data=[1, 2, 3]) == T(data=[1, 2, 3])
                        )

                        tref = T()
                        tref.b = b
                        tref.delta = delta
                        tref.t = tr
                        tref.evt = evt
                        tref.sta = sta
                        tref.picks = picks
                        tref.meta = meta
                        @test T(; b, delta, data=tr, evt, sta, picks, meta) == tref
                    end

                    @testset "Parameterised constructor" begin
                        @test T{Float32}(b=1, delta=0.5, data=[1,2,3]) ==
                            T{Float32}(1, 0.5, [1,2,3])
                        @test T{Float32,Vector{Float16}}(b=1, delta=0.5, data=[1,2,3]) ==
                            T{Float32,Vector{Float16}}(1, 0.5, [1,2,3])
                    end
                end

                @testset "Parameterised constructor" begin
                    @test typeof(Trace{Float16}()) == Trace{Float16,Vector{Float16},Seis.Geographic{Float16}}
                    @test typeof(Trace{Float16,Vector{Float32}}()) ==
                        Trace{Float16,Vector{Float32},Seis.Geographic{Float16}}
                    @test typeof(Trace{Float16,Vector{Float32},Seis.Cartesian{Float16}}()) ==
                        Trace{Float16,Vector{Float32},Seis.Cartesian{Float16}}
                end
            end
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
        let t = Trace{Float32,Vector{Float32}}(rand(), rand(), rand(3))
            @test t isa Trace{Float32,Vector{Float32},Seis.Geographic{Float32}}
        end
        @test Trace{Float32}(0, 1, 100) isa Trace{Float32,Vector{Float32},Seis.Geographic{Float32}}

        @test_throws TypeError Trace{Int,Vector{Int}}(1, 1, rand(1:3, 3))
    end

    @testset "FourierTrace" begin
        @testset "Explicit constructor" begin
            @testset "T = $T" for T in (Float32, Float64)
                @testset "V = Vector{Complex{$V}}" for V in (Float32, Float64)
                    @testset "Geom = $P" for P in (Seis.Geographic{T}, Seis.Cartesian{T})
                        b = rand(T)
                        delta = rand(T)
                        n = rand(1:100)
                        tr = rand(Complex{V}, n)
                        f = FourierTrace{T,Vector{Complex{V}},P}(;
                            b=b, delta=delta, data=tr)
                        @test f isa FourierTrace{T,Vector{Complex{V}},P}
                        @test nfrequencies(f) == n
                        @test nsamples(f) == 2*nfrequencies(f) - 1
                        @test all(tr .== trace(f))
                        @test f.b == b
                        @test f.delta == delta
                        @test f.evt isa Event{T,P}
                        @test f.sta isa Station{T,P}
                    end
                end
            end

            @testset "Errors" begin
                @test_throws UndefKeywordError FourierTrace()
                @test_throws UndefKeywordError FourierTrace(; b=1)
                @test_throws UndefKeywordError FourierTrace(; b=1, delta=0.1)
                @test_throws ArgumentError FourierTrace(;
                    b=1, delta=-0.1, data=rand(Complex{Float64}, 10))
                @test_throws ArgumentError FourierTrace(;
                    b=1, delta=0.1, data=rand(Complex{Float64}, 10), nsamples=-1)
                @test_throws MethodError FourierTrace(; b=1, delta=0.1, data=1)
            end

            @testset "Default parameters" begin
                b = 0
                delta = 1
                data = rand(Complex{Float32}, 10)
                @test (FourierTrace(; b=b, delta=delta, data=data)
                    isa FourierTrace{Float64,Vector{Complex{Float64}},Seis.Geographic{Float64}})

                @testset "T = $T" for T in (Float32, Float64)
                    @test (FourierTrace{T}(; b=0, delta=1, data=data) isa
                        FourierTrace{T,Vector{Complex{T}},Seis.Geographic{T}})

                    @testset "V = Vector{Complex{$V}}" for V in (Float32, Float64)
                        @test (FourierTrace{T,Vector{Complex{V}}}(; b, delta, data)
                            isa FourierTrace{T,Vector{Complex{V}},Seis.Geographic{T}})
                    end
                end
            end
        end

        @testset "fft constructor" begin
            b = rand()
            delta = rand()
            n = rand(1:100)
            data = rand(n)
            t = Trace(b, delta, data)
            f = fft(t)
            fdelta = 1/(n*delta)
            @test f == FourierTrace(; b=b, nsamples=n, delta=fdelta, data=FFTW.rfft(data))
        end

        @testset "ifft constructor" begin
            @testset "Original number of points" begin
                f = fft(Trace(0, 1, rand(3)))
                @test nsamples(ifft(f)) == 3
                @test nsamples(ifft(f, 3)) == 3
                @test nsamples(ifft(f, 2)) == 3
            end

            @testset "Modified trace length" begin
                f = fft(Trace(0, 1, rand(12)))
                foreach(_ -> pop!(trace(f)), 1:3)
                @test nsamples(ifft(f)) == 6
                @test nsamples(ifft(f, 6)) == 6
                @test nsamples(ifft(f, 7)) == 7
            end
        end

        @testset "Accessors" begin
            @testset "$T" for T in (Float32, Float64)
                @testset "Vector{$V}" for V in (Float32, Float64)
                    b = rand(T)
                    delta = rand(T)
                    n = rand(5:100)
                    data = rand(V, n)
                    t = Trace{T,Vector{V}}(b, delta, data)
                    f = fft(t)
                    @test f isa FourierTrace{T,Vector{Complex{V}},Seis.Geographic{T}}
                    @test nsamples(f) == n
                    @test trace(f) == FFTW.rfft(data)
                    @test nfrequencies(f) == length(trace(f)) == (n÷2) + 1
                    @test frequencies(f) == (0:(length(trace(f)) - 1)).*f.delta
                    @test starttime(f) == b
                    @test eltype(f) == Complex{V}
                end
            end

            @testset "nsamples" begin
                @test nsamples(fft(Trace(0, 1, rand(3)))) == 3
                @test nsamples(fft(Trace(0, 1, rand(3))), even=true) == 2
                @test nsamples(fft(Trace(0, 1, rand(3))), even=false) == 3

                @test nsamples(FourierTrace(
                    b=0, delta=1, data=complex.(rand(3)))) == 5
                @test nsamples(FourierTrace(
                    b=0, delta=1, data=complex.(rand(3))), even=true) == 4
                @test nsamples(FourierTrace(
                    b=0, delta=1, data=complex.(rand(3))), even=false) == 5
            end
        end

        @testset "Type" begin
            let b = 1, delta = 1, data = complex.(rand(1:3, 100), 0)
                f = FourierTrace(; b, delta, data)
                @test f isa FourierTrace
                @test typeof(f) <: Seis.AbstractFourierTrace
                @test typeof(f) <: FourierTrace
                @test typeof(f) == FourierTrace{Float64,Vector{Complex{Float64}},Seis.Geographic{Float64}}
            end
        end

        @testset "Arrays" begin
            @testset "Broadcasting as scalar" begin
                @test nsamples.(Trace(0, 1, rand(2))) == 2
            end

            @testset "getproperty/setproperty!" begin
                b = rand()
                t = FourierTrace(; b=b, delta=rand(), data=complex.(rand(3)))
                ts = [deepcopy(t) for _ in 1:3]
                @test ts isa Vector{typeof(t)}
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
                    kwargs = Geom == Seis.Cartesian ? (x=1, y=2, z=3) :
                        ES == Station ? (lon=1, lat=2, elev=3) : (lon=1, lat=2, dep=-0.003)
                    if ES == Station
                        kwargs = (net="A", sta="B", loc="C", cha="D", elev=-1,
                                  azi=rand(), inc=rand(), meta=Dict(:x=>"X"), kwargs...)
                    else
                        kwargs = (id="XYZ", time=now(), meta=Dict(:y=>"Y"), kwargs...)
                    end
                    @test convert(outtype, intype()) isa outtype
                    @test convert(outtype, intype(; kwargs...)) isa outtype
                    @test getfield(convert(outtype, intype(; kwargs...)), :pos) == Geom{Tout}(1, 2, 3)
                    # Identity 'conversions' give the same object
                    if intype == outtype
                        in = intype(; kwargs...)
                        @test convert(outtype, in) === in
                    end
                end

                @testset "$Tr $Geom" for Tr in (Trace, FourierTrace)
                    Vs = if Tr == FourierTrace
                        [Vector{Complex{TT}} for TT in Ts]
                    else
                        [Vector{TT} for TT in Ts]
                    end
                    @testset "Data type $Vin -> $Vout" for Vin in Vs, Vout in Vs
                        intype = Tr{Tin, Vin, Geom{Tin}}
                        outtype = Tr{Tout, Vout, Geom{Tout}}
                        vtype = eltype(Vout)
                        data = vtype[1, 2, 3]
                        tin = if Tr == FourierTrace
                            intype(; b=0, delta=1, data=data)
                        else
                            intype(0, 1, data)
                        end
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

        @testset "FourierTrace" begin
            fs = [fft(Trace(rand(), 1, [1, 2, 3])) for _ in 1:3]
            fv = view(fs, 2:3)
            @test fv.b == fs[2:3].b
            @test fv[1] == fs[2]
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
        @testset "FourierTrace" begin
            f = fft.([Trace(0, 1, [1,2,3]) for _ in 1:3])
            push!(trace(f[3]), 0)
            test_hash_isequal(f...)
            pop!(trace(f[3]))
            f[3].b = 1
            test_hash_isequal(f...)
            f[3].b = f[1].b
            f[3].sta.lon = 1
            test_hash_isequal(f...)
            f[3].sta.lon = missing
            f[3].evt.time = now()
            test_hash_isequal(f...)
            f[3].evt.time = missing
            f[3].meta.a = 1
            test_hash_isequal(f...)
        end
    end
end
