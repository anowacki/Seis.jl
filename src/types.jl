# Types used by Seis

"Default floating point type used for traces if none specified as type parameters"
const DEFAULT_FLOAT = Float64
"Default string type used for traces if none specified as type parameters"
const DEFAULT_STRING = String

"""
    SeisDict

Wrapper around `Base.Dict` which allows one to get or set values using
the {get|set}property[!] syntax, i.e., `.`-access like `dict.key`.
`SeisDict`s also differ in that `missing` is returned instead of throwing a
`KeyError` when accessing a nonexistent key.  A key is removed if it is set to
`missing`.

```
julia> dict = Seis.SeisDict(:a=>1)
Dict{Any,Any} with 1 entry:

julia> dict.a
1

julia> dict.b
missing
```

Arrays of `SeisDict`s also define `.`-access and setting via broadcasting, so one
may do:

```
julia> d1 = Seis.SeisDict(:a=>1, :b=>2)
Dict{Any,Any} with 2 entries:
  :a => 1
  :b => 2

julia> d2 = deepcopy(d1); d.a = 2;

julia> d = [d1, d2];

julia> d.a
2-element Array{Int64,1}:
 1
 1
```
"""
struct SeisDict{K,V} <: Base.AbstractDict{K,V}
    dict::Dict{K,V}
end

SeisDict{K,V}(args...) where {K,V} = SeisDict{K,V}(Dict{K,V}(args...))
SeisDict{K}(args...) where {K} = SeisDict{K,Any}(Dict{K,Any}(args...))
SeisDict(args...) = SeisDict{Any,Any}(Dict{Any,Any}(args...))
SeisDict(dict::Dict{K,V}) where {K,V} = SeisDict{K,V}(dict)

# Have to use getfield here because we define getproperty as getindex
Base.getindex(sd::SeisDict{K,V}, key) where {K,V} = get(getfield(sd, :dict), key, missing)
Base.setindex!(sd::SeisDict{K,V}, ::Missing, key) where {K,V} = delete!(getfield(sd, :dict), key)
Base.setindex!(sd::SeisDict{K,V}, v, i...) where {K,V} = setindex!(getfield(sd, :dict), v, i...)
Base.Dict(sd::SeisDict) = getfield(sd, :dict)
Base.summary(sd::SeisDict) = summary(getfield(sd, :dict))
Base.iterate(sd::SeisDict, args...) = iterate(getfield(sd, :dict), args...)
Base.length(sd::SeisDict) = length(getfield(sd, :dict))
Base.get(sd::SeisDict, args...) = get(getfield(sd, :dict), args...)
Base.delete!(sd::SeisDict, args...) = (delete!(getfield(sd, :dict), args...); sd)

Base.getproperty(sd::SeisDict, key::Symbol) = getindex(sd, key)
Base.setproperty!(sd::SeisDict, key::Symbol, val) = setindex!(sd, val, key)
Base.propertynames(sd::SeisDict, private=false) = collect(keys(sd))

Base.getproperty(sd::AbstractArray{<:SeisDict}, key::Symbol) = getindex.(sd, key)
Base.setproperty!(sd::AbstractArray{<:SeisDict}, key::Symbol, val) = setindex!.(sd, val, key)

"""
    Event

Type containing information about a seismic event.  Fields `lon` and `lat` are
epicentral location in °; `dep` is depth below the reference (e.g., sea level) in km.
`time` is a `Dates.DateTime` giving the event origin date and time, while `id`
is a string holding the event identifier.
`meta` is a `Dict` holding any extra information about the event.

Missing information is allowed and stored as `missing`.
"""
mutable struct Event{T<:AbstractFloat,S<:AbstractString}
    lon::Union{T,Missing}
    lat::Union{T,Missing}
    dep::Union{T,Missing}
    time::Union{DateTime,Missing}
    id::Union{S,Missing}
    meta::SeisDict{Symbol,Any}
end
Event{T,S}() where {T,S} = Event{T,S}(missing, missing, missing, missing, missing, Dict())
Event(lon=missing, lat=missing, dep=missing, time=missing, kind=missing, meta=Dict()) =
    Event{DEFAULT_FLOAT,DEFAULT_STRING}(lon, lat, dep, time, kind, meta)

const EVENT_FIELDS = fieldnames(Event)
Base.getproperty(e::AbstractArray{<:Event}, f::Symbol) = getfield.(e, f)
Base.setproperty!(e::AbstractArray{<:Event}, f::Symbol, val) =
    setfield!.(e, f, convert.(typeof.(getfield.(e, f)), val))
Base.propertynames(e::AbstractArray{<:Event}, private=false) = fieldnames(eltype(e))
Base.:(==)(e1::Event, e2::Event) =
    all(x -> isequal(x[1], x[2]), (getfield.((e1, e2), f) for f in EVENT_FIELDS))

"""
    Station

Struct containing information about a seismic station.  Fields `net`, `sta`, `loc`
and `cha` are the station, network, channel and location codes respectively, whilst
`lon` and `lat` are the location in °.  Set station depth `dep` and elevation `elev`
in m relative to the reference level.  The azimuth `azi` and inclination `inc` of
the channel in ° are respectively measured from north to east, and downward from the
vertical.  (E.g., a "BHN" channel typically will have a `azi == 0` and `inc == 90`.)

`meta` is a `Dict` holding any extra
information about the station or channel.

Missing information is allowed and stored as `missing`.
"""
mutable struct Station{T<:AbstractFloat,S<:AbstractString}
    net::Union{S,Missing}
    sta::Union{S,Missing}
    loc::Union{S,Missing}
    cha::Union{S,Missing}
    lon::Union{T,Missing}
    lat::Union{T,Missing}
    dep::Union{T,Missing}
    elev::Union{T,Missing}
    azi::Union{T,Missing}
    inc::Union{T,Missing}
    meta::SeisDict{Symbol,Any}
end
Station{T,S}() where {T,S} = Station{T,S}(missing, missing, missing, missing, missing,
                                          missing, missing, missing, missing, missing,
                                          Dict())
Station(net=missing, sta=missing, loc=missing, cha=missing, lon=missing, lat=missing,
        dep=missing, elev=missing, azi=missing, inc=missing, meta=Dict()) =
    Station{DEFAULT_FLOAT,DEFAULT_STRING}(net, sta, loc, cha, lon, lat, dep, elev,
                                          azi, inc, meta)

const STATION_FIELDS = fieldnames(Station)
Base.getproperty(s::AbstractArray{<:Station}, f::Symbol) = getfield.(s, f)
Base.setproperty!(s::AbstractArray{<:Station}, f::Symbol, val) =
    setfield!.(s, f, convert.(typeof.(getfield.(s, f)), val))
Base.propertynames(s::AbstractArray{<:Station}, private=false) = fieldnames(eltype(s))
Base.:(==)(s1::Station, s2::Station) =
    all(x -> isequal(x[1], x[2]), (getfield.((s1, s2), f) for f in STATION_FIELDS))

"""
    Pick{T,S}

`NamedTuple` defining a named arrival time which can be associated with a `Trace`.

If `p` is a `Pick`, then `p.time` gives the arrival time, and `p.name` gives its
description.

Arrays of `Pick`s, which are stored in `Traces`, can be accessed using
`getproperty` (`.`-access) and this returns arrays of times or names.
"""
const Pick{T,S} = NamedTuple{(:time, :name), Tuple{T,Union{Missing,S}}} where {T,S}
Base.getproperty(p::AbstractVector{<:Pick}, field::Symbol) = getproperty.(p, field)

"""
    AbstractTrace

Abstract type from which you should subtype if creating new types of traces.
"""
abstract type AbstractTrace end

"""
    Trace

Evenly-sampled time series recorded at a single seismic station.  The start time
of the trace, in s, is in the `b` field, whilst the sampling interval, in s, is `delta`.
The trace itself is accessed using the `trace` method, like `trace(t)`.

All `Trace`s are relative to the event time `evt.time` if it is defined, regardless
of what the event is.  For example, `evt.time` could be the origin time of an
earthquake, or a picked arrival time.

The trace then contains information about an associated `Event` in `evt` and the
`Station` in `sta`.  `picks` holds a vector which contains pairs of pick times
relative to the origin and names (which can be `missing`).  Access picks with the
`picks` method, and add picks with `add_pick!`.

The `meta` `Dict` holds any other information about the trace.

If the event `time` is set, then the trace beginning time `b` is relative to this.
"""
mutable struct Trace{T<:AbstractFloat,V<:AbstractVector{T},S<:AbstractString} <: AbstractTrace
    b::T
    delta::T
    t::V
    evt::Event{T,S}
    sta::Station{T,S}
    picks::Vector{Pick{T,S}}
    meta::SeisDict{Symbol,Any}
    function Trace{T,S,V}(b, delta, t, evt, sta, picks, meta) where {T,S,V}
        delta > 0 || throw(ArgumentError("delta cannot be <= 0"))
        new(b, delta, t, evt, sta, picks, meta)
    end
end

Trace{T,V,S}(b, delta, t) where {T,V,S} = Trace{T,V,S}(b, delta, t, Event{T,S}(),
                                                       Station{T,S}(), [], Dict())

function Trace(b::T1, delta::T2, t::AbstractVector{T3}) where {T1,T2,T3}
    T = float(promote_type(T1, T2, T3))
    Trace{T,Vector{T},DEFAULT_STRING}(b, delta, t, Event{T,DEFAULT_STRING}(),
                                      Station{T,DEFAULT_STRING}(), [], Dict())
end

const TRACE_FIELDS = fieldnames(Trace)

Base.getproperty(t::AbstractArray{<:Trace}, f::Symbol) = getfield.(t, f)
Base.setproperty!(t::AbstractArray{<:Trace}, f::Symbol, val) =
    setfield!.(t, f, convert.(typeof.(getfield.(t, f)), val))
Base.propertynames(t::AbstractArray{<:Trace}, private=false) = fieldnames(eltype(t))
Base.:(==)(t1::Trace, t2::Trace) =
    all(x -> isequal(x[1], x[2]), (getfield.((t1, t2), f) for f in TRACE_FIELDS))
