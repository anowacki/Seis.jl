# Types used by Seis

"Default floating point type used for traces if none specified as type parameters"
const DEFAULT_FLOAT = Float64

"""
    SeisDict

Wrapper around `Base.Dict` which allows one to get or set values using
the {get|set}property[!] syntax, i.e., `.`-access like `dict.key`.
`SeisDict`s also differ in that `missing` is returned instead of throwing a
`KeyError` when accessing a nonexistent key.  A key is removed if it is set to
`missing`.

# Examples

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

julia> d2 = deepcopy(d1); d2.a = 2;

julia> d = [d1, d2];

julia> d.a
2-element Array{Int64,1}:
 1
 2
```
"""
struct SeisDict{K,V} <: Base.AbstractDict{K,V}
    dict::Dict{K,V}
end

SeisDict{K,V}(args...) where {K,V} = SeisDict{K,V}(Dict{K,V}(args...))
SeisDict{K}(args...) where {K} = SeisDict{K,Any}(Dict{K,Any}(args...))
SeisDict(args...) = SeisDict{Any,Any}(Dict{Any,Any}(args...))

# Have to use getfield here because we define getproperty as getindex
Base.getindex(sd::SeisDict{K,V}, key) where {K,V} = get(getfield(sd, :dict), key, missing)
Base.setindex!(sd::SeisDict{K,V}, ::Missing, key) where {K,V} = delete!(getfield(sd, :dict), key)
Base.setindex!(sd::SeisDict{K,V}, v, i...) where {K,V} = setindex!(getfield(sd, :dict), v, i...)
Base.Dict(sd::SeisDict) = getfield(sd, :dict)
Base.empty!(sd::SeisDict) = empty!(getfield(sd, :dict))
Base.summary(sd::SeisDict) = summary(getfield(sd, :dict))
Base.iterate(sd::SeisDict, args...) = iterate(getfield(sd, :dict), args...)
Base.length(sd::SeisDict) = length(getfield(sd, :dict))
Base.get(sd::SeisDict, args...) = get(getfield(sd, :dict), args...)
Base.delete!(sd::SeisDict, args...) = (delete!(getfield(sd, :dict), args...); sd)

Base.getproperty(sd::SeisDict, key::Symbol) = getindex(sd, key)
Base.setproperty!(sd::SeisDict{K,V}, key::Symbol, val) where {K,V} = setindex!(sd, V(val), key)
Base.setproperty!(sd::SeisDict{K,Any}, key::Symbol, val) where K = setindex!(sd, val, key)
Base.propertynames(sd::SeisDict, private=false) = collect(keys(sd))

Base.getproperty(sd::AbstractArray{<:SeisDict}, key::Symbol) = getindex.(sd, key)
Base.setproperty!(sd::AbstractArray{<:SeisDict}, key::Symbol, val) = setindex!.(sd, val, key)

"""
    Position

Abstract type of which other position types are subtypes.
A `Position` specifies where in space an object is located.

The coordinates of a `Position` can be accessed via `getproperty`
(e.g., `p.x`) or index (e.g., `p[1]`):

# Examples
```
julia> p = Seis.Cartesian(x=1, y=2, z=3)
Seis.Cartesian{Float64}(1.0, 2.0, 3.0)

julia> p.x === p[1]
true
```
"""
abstract type Position{T} end

Base.getproperty(p::AbstractArray{<:Position}, f::Symbol) = getfield.(p, f)
Base.getindex(p::Position, i::Int) = getfield(p, i)
Base.:(==)(p1::P, p2::P) where {P<:Position} =
    all(isequal(getfield(p1, f), getfield(p2, f)) for f in fieldnames(P))

"Return the parameter `T` for a subtype of `Position{T}`"
@inline _eltype(::Type{<:Position{T}}) where T = T

"""
    Geographic{T} <: Position

A geographic position in spherical coordinates.  Accessible fields are:

- `lon`: Longitude (°)
- `lat`: Latitude (°)
- `dep`: Depth below the reference level (e.g., ellipsoid) (km)

It is recommended that for `Station`s, the `dep` field describes the
depth of the sensor relative to sea level in km, so a station at
150 m elevation has a depth of –0.15 km.  Information about sensor
burial depth should be held in the `Event`'s `meta` field.
"""
mutable struct Geographic{T<:AbstractFloat} <: Position{T}
    lon::Union{Missing,T}
    lat::Union{Missing,T}
    dep::Union{Missing,T}
end
Geographic{T}(; lon=missing, lat=missing, dep=missing) where T =
    Geographic{T}(lon, lat, dep)
Geographic(args...; kwargs...) = Geographic{DEFAULT_FLOAT}(args...; kwargs...)

"""
    Cartesian{T} <: Position{T}

A position in Cartesian coordinates.  Accessible fields are:

- `x`: X coordinate (m)
- `y`: Y coordinate (m)
- `z`: Z coordinate (m)

Seis.jl's convention is that `x` and `y` are horizontal (usually with `x`
being the local easting and `y` the local northing),
and `z` is upwards (giving a right-handed system).  The units are m.
"""
mutable struct Cartesian{T<:AbstractFloat} <: Position{T}
    x::Union{Missing,T}
    y::Union{Missing,T}
    z::Union{Missing,T}
end
Cartesian{T}(; x=missing, y=missing, z=missing) where T =
    Cartesian{T}(x, y, z)
Cartesian(args...; kwargs...) = Cartesian{DEFAULT_FLOAT}(args...; kwargs...)

"""
    Event

Type containing information about a seismic event.  Fields `lon` and `lat` are
epicentral location in °; `dep` is depth below the reference (e.g., sea level) in km.
`time` is a `Dates.DateTime` giving the event origin date and time, while `id`
is a string holding the event identifier.
`meta` is a `Dict` holding any extra information about the event.

Missing information is allowed and stored as `missing`.
"""
mutable struct Event{T<:AbstractFloat, P<:Position}
    pos::P
    time::Union{DateTime,Missing}
    id::Union{String,Missing}
    meta::SeisDict{Symbol,Any}
end

"""
    Event(; kwargs...)

Construct a type representing a seismic event (which can simply be a
start time for the trace with no extra information).

## Keyword arguments
All arguments default to `missing`.
- `time`: `DateTime` of the reference time
- `id`: string giving the event identifier
- `meta`: `Dict` holding any extra information

By default, the following may be passed:
- `lon`: event longitude in °
- `lat`: event latitude in °
- `dep`: event depth in km

---
    Event{T,P}(; kwargs...)
    Event{T}(; kwargs...)

Construct an `Event` with specific floating point type `T` and
`Position` type `P`.  If not given, `T` defaults to `$DEFAULT_FLOAT`
and `P` defaults to `Seis.Geographic{$DEFAULT_FLOAT}`.

See also [`CartEvent`](@ref).
"""
Event{T,P}(; lon=missing, lat=missing, dep=missing, time=missing,
        id=missing, meta=SeisDict{Symbol,Any}()) where {T, P<:Geographic} =
    Event{T,P}(P(lon, lat, dep), time, id, meta)

Event{T,P}(; x=missing, y=missing, z=missing, time=missing,
        id=missing, meta=SeisDict{Symbol,Any}()) where {T, P<:Cartesian} =
    Event{T,P}(P(x, y, z), time, id, meta)

Event{T}(; kwargs...) where T = Event{T ,Geographic{T}}(; kwargs...)
Event(; kwargs...) = Event{DEFAULT_FLOAT, Geographic{DEFAULT_FLOAT}}(; kwargs...)

"""
    CartEvent{T}

Alias for `Event{T, Cartesian{T}} where T`, representing an [`Event`](@ref)
with Cartesian coordinates.

This type is useful for dispatch, allowing one to write methods which
are only applicable when a `Event` has Cartesian coordinates.

# Example
```
julia> using Geodesy

julia> Geodesy.LLA(evt::CartEvent) = LLA(evt.x, evt.y, evt.z)
```
"""
const CartEvent{T} = Event{T,Cartesian{T}}

"""
    CartEvent{T}(; kwargs...)
    CartEvent(; kwargs...)

Construct a `CartEvent`.  See [`Station`](@ref) for details on keyword
arguments, noting that position must be given via a combination of the
keyword arguments `x`, `y` and `z` (not `lon`, `lat` or `dep`).
"""
CartEvent(; kwargs...) = Event{DEFAULT_FLOAT, Cartesian{DEFAULT_FLOAT}}(; kwargs...)

"""
    GeogEvent{T} where T

Alias for `Event{T, Geographic{T}} where T`, representing an [`Event`](@ref)
with geographic coordinates.

This type is useful for dispatch, allowing one to write methods which
are only applicable when a `Event` has geographic coordinates.

Note that `Event`s are geographic by default.  Use [`Event()`](@ref Event)
to construct a geographic `Event`.

# Example
```
julia> using Geodesy

julia> Geodesy.LLA(evt::GeogEvent) = LLA(evt.lat, evt.lon, evt.elev)
```
"""
const GeogEvent{T} = Event{T, Geographic{T}}

const EVENT_FIELDS = fieldnames(Event)

Base.getproperty(e::AbstractArray{<:Event}, f::Symbol) = getproperty.(e, f)
Base.setproperty!(e::AbstractArray{<:Event}, f::Symbol, val) =
    setproperty!.(e, f, val)
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
mutable struct Station{T<:AbstractFloat, P<:Position{T}}
    net::Union{String,Missing}
    sta::Union{String,Missing}
    loc::Union{String,Missing}
    cha::Union{String,Missing}
    pos::P
    elev::Union{T,Missing}
    azi::Union{T,Missing}
    inc::Union{T,Missing}
    meta::SeisDict{Symbol,Any}
end

"""
    Station(; kwargs...)

Construct a type representing a seismic station.

## Keyword arguments
All values default to `missing`.
- `net`: network code
- `sta`: station code
- `loc`: location code
- `cha`: channel code
- `elev`: local station elevation above the ground in m
- `azi`: component azimuth, in ° east of local north
- `inc`: component inclination, in ° down from upwards
- `meta`: `Dict` holding any extra information

By default, the following may be passed:
- `lon`: longitude in °
- `lat`: latitude in °
- `dep`: depth in km

If `P <: Seis.Cartesian` (see [`Station{T,S,P}`](@ref) or
[`CartStation`](@ref)), then the following may be passed:
- `x`: X coordinate in m
- `y`: Y coordinate in m
- `z`: Z coordinate in m

---
    Station{T,P}(; kwargs...)
    Station{T}(; kwargs...)

Construct a `Station` with specific floating point type `T` and
`Position` type `P`.  If not given, `T` defaults to `$DEFAULT_FLOAT`
and `P` defaults to `Seis.Geographic{$DEFAULT_FLOAT}`.
"""
Station{T,P}(; net=missing, sta=missing, loc=missing,
        cha=missing, lon=missing, lat=missing, dep=missing, elev=missing,
        azi=missing, inc=missing, meta=Dict()) where {T, P<:Geographic} =
    Station{T, P}(net, sta, loc, cha, P(lon, lat, dep), elev, azi, inc, meta)

Station{T,P}(; net=missing, sta=missing, loc=missing,
        cha=missing, x=missing, y=missing, z=missing, elev=missing,
        azi=missing, inc=missing, meta=Dict()) where {T, P<:Cartesian} =
    Station{T,P}(net, sta, loc, cha, P(x, y, z), elev, azi, inc, meta)

Station{T}(; kwargs...) where T = Station{T, Geographic{T}}(; kwargs...)
Station(; kwargs...) = Station{DEFAULT_FLOAT, Geographic{DEFAULT_FLOAT}}(; kwargs...)

"""
    CartStation{T} where T

Alias for `Station{T, Cartesian{T}} where T`, representing a
[`Station`](@ref) with Cartesian coordinates.

This type is useful for dispatch, allowing one to write methods which
are only applicable when a `Station` has Cartesian coordinates.

# Example
Create a function which obtains the east-north-up coordinates of
a `Station`
```
julia> using Geodesy

julia> Geodesy.ENU(sta::CartStation) = ENU(sta.x, sta.y, sta.z)

julia> ENU(CartStation(x=1, y=2, z=3))
```
"""
const CartStation{T} = Station{T, Cartesian{T}}

"""
    CartStation{T,S}(; kwargs...)
    CartStation{T}(; kwargs...)
    CartStation(; kwargs...)

Construct a `CartStation`.  See [`Station`](@ref) for details on keyword
arguments, noting that position must be given via a combination of the
keyword arguments `x`, `y` and `z` (not `lon`, `lat` or `dep`).
"""
CartStation(; kwargs...) = CartStation{DEFAULT_FLOAT}(; kwargs...)

Base.propertynames(e::Union{Event{T,P}, Station{T,P}}) where {T,P} =
    (fieldnames(P)..., fieldnames(typeof(e))...)

"""
    GeogStation{T} where T

Alias for `Station{T, Geographic{T}} where T`, representing a [`Station`](@ref)
with geographic coordinates.

This type is useful for dispatch, allowing one to write methods which
are only applicable when a `Station` has geographic coordinates.

Note that `Station`s are geographic by default.  Use [`Station()`](@ref Station)
to construct a geographic `Station`.

# Example
```
julia> using Geodesy

julia> Geodesy.LLA(sta::GeogStation) = LLA(sta.lat, sta.lon, sta.elev)
"""
const GeogStation{T} = Station{T, Geographic{T}}

const STATION_FIELDS = fieldnames(Station)
function Base.getproperty(e::Union{Station,Event}, p::Symbol)
    if p === :lon || p === :lat || p === :dep || p === :x || p === :y || p === :z
        getfield(e.pos, p)
    else
        getfield(e, p)
    end
end

function Base.setproperty!(e::Union{Event{T},Station{T}}, p::Symbol, v) where T
    if p === :lon || p === :lat || p === :dep || p === :x || p === :y || p === :z
        setfield!(e.pos, p, convert(Union{Missing,T}, v))
    else
        # setfield!(e, p, v)
        setfield!(e, p, convert(fieldtype(typeof(e), p), v))
    end
end

Base.getproperty(s::AbstractArray{<:Station}, f::Symbol) = getproperty.(s, f)
Base.setproperty!(s::AbstractArray{<:Station}, f::Symbol, val) =
    setproperty!.(s, f, val)
Base.propertynames(s::AbstractArray{<:Station}, private=false) = propertynames(first(s))
Base.:(==)(s1::Station, s2::Station) =
    all(x -> isequal(x[1], x[2]), (getfield.((s1, s2), f) for f in STATION_FIELDS))

"""
    Pick{T}

`NamedTuple` defining a named arrival time which can be associated with a `Trace`.

If `p` is a `Pick`, then `p.time` gives the arrival time, and `p.name` gives its
description.

Arrays of `Pick`s, which are stored in `Traces`, can be accessed using
`getproperty` (`.`-access) and this returns arrays of times or names.

# Example
```
julia> t = Trace(0, 1, 1)

julia> t.picks.P = 1.2

julia> t.picks.S = 2.1, "Sn"

julia> t.picks.P.time

julia> t.picks
"""
const Pick{T} = NamedTuple{(:time, :name), Tuple{T,Union{Missing,String}}} where T
Base.getproperty(p::AbstractVector{<:Pick}, field::Symbol) = getproperty.(p, field)
# Allow conversion from a single time
Pick{T}(time::Real) where T = Seis.Pick{T}((time=time, name=missing))

"""
    AbstractTrace

Abstract type from which you should subtype if creating new types of traces.
"""
abstract type AbstractTrace end

"""
    Trace

Evenly-sampled time series recorded at a single seismic station.  The start time
of the trace, in s, is in the `b` property, whilst the sampling interval, in s,
is `delta`.
The trace itself is accessed using the `trace` method, like `trace(t)`.

All `Trace`s are relative to the event time `evt.time` if it is defined, regardless
of what the event is.  For example, `evt.time` could be the origin time of an
earthquake, or a picked arrival time.

The trace then contains information about an associated `Event` in `evt` and the
`Station` in `sta`.  `picks` holds a dictionary which contains pairs of pick times
relative to the origin and names (which can be `missing`).  Access picks with the
`picks` method, and add picks with `add_pick!`.

The `meta` `Dict` holds any other information about the trace.

If the event `time` is set, then the trace beginning time `b` is relative to this.

Find the trace start time relative to the origin time using [`starttime`](@ref).
The absolute start time and date, if an origin time is set, is given by
[`startdate`](@ref).
"""
mutable struct Trace{T<:AbstractFloat,V<:AbstractVector{<:AbstractFloat},P<:Position{T}} <: AbstractTrace
    b::T
    delta::T
    t::V
    evt::Event{T,P}
    sta::Station{T,P}
    picks::SeisDict{Union{Symbol,Int},Pick{T}}
    meta::SeisDict{Symbol,Any}
    function Trace{T,V,P}(b, delta, t, evt, sta, picks, meta) where {T,V,P}
        delta > 0 || throw(ArgumentError("delta cannot be <= 0"))
        new(b, delta, t, evt, sta, picks, meta)
    end
end

"""
    Trace(b, delta, t::AbstractVector) -> trace::Trace{$DEFAULT_FLOAT,Vector{$DEFAULT_FLOAT},Seis.Geographic{$DEFAULT_FLOAT}}

Create a `Trace` with starting time `b` s, sampling interval `delta` s and
an `AbstractVector` `t` of values for the trace.  The default precision for
the type is `$DEFAULT_FLOAT`.  By default, the trace's event and station
are in geographic coordinates.

    Trace(b, delta, n::Integer) -> trace::Trace{$DEFAULT_FLOAT,Vector{$DEFAULT_FLOAT},Seis.Geographic{$DEFAULT_FLOAT}}

Create a new `Trace` with uninitialised data of length `n` samples.

    Trace{T,V,P}(args...) -> trace::Trace{T,V,P}
    Trace{T,V}(args...) -> trace::Trace{T,V,Seis.Geographic{T}}
    Trace{T}(args...) -> trace::Trace{T,Vector{T},Seis.Geographic{T}}

Construct traces with non-default precision and string types.  In the third
form, the vector type defaults to `Vector{$DEFAULT_FLOAT}`, whilst in
both the second and third forms, `P` defaults to `Seis.Geographic{T}`.
"""
Trace{T,V,P}(b, delta, t::AbstractVector) where {T,V,P} =
    Trace{T,V,P}(b, delta, t, Event{T,P}(), Station{T,P}(), Dict(), Dict())
Trace{T,V,P}(b, delta, n::Integer) where {T,V,P} = Trace{T,V,P}(b, delta, Vector{T}(undef, n))
Trace{T,V}(args...) where {T,V} = Trace{T,V,Geographic{T}}(args...)
Trace{T}(args...) where T = Trace{T,Vector{T},Geographic{T}}(args...)
Trace(b, delta, t_or_n) = Trace{DEFAULT_FLOAT,Vector{DEFAULT_FLOAT},Geographic{DEFAULT_FLOAT}}(b, delta, t_or_n)

"""
    CartTrace

Alias for `Trace` where `Event` and `Station` coordinates are
`Seis.Cartesian` rather than `Seis.Geographic`
"""
const CartTrace{T,V} = Trace{T, V, Cartesian{T}}

"""
    CartTrace

Create a new `Trace` whose `Event` and `Station` are in Cartesian coordinates.

See [`Trace`](@ref) for more details of the different construction methods.

---
    CartTrace{T,V}(b, delta, t::AbstractVector) -> trace
    CartTrace{T}(args...)

Construct a `Trace` with Cartesian coordinates for the `Event` and `Station`
with non-default precision and string types.  See [`Trace{T,V}`](@ref)
for details.
"""
CartTrace{T}(args...) where T = Trace{T, Vector{T}, Cartesian{T}}(args...)
CartTrace(args...) = Trace{DEFAULT_FLOAT, Vector{DEFAULT_FLOAT},
          Cartesian{DEFAULT_FLOAT}}(args...)

const TRACE_FIELDS = fieldnames(Trace)

Base.getproperty(t::AbstractArray{<:Trace}, f::Symbol) = getfield.(t, f)
Base.setproperty!(t::AbstractArray{<:Trace}, f::Symbol, val) =
    for tt in t setproperty!(tt, f, val) end
function Base.setproperty!(t::AbstractArray{<:Trace}, f::Symbol, val::AbstractArray)
    length(t) == length(val) || throw(DimensionMismatch())
    for (tt, vv) in zip(t, val)
        setproperty!(tt, f, vv)
    end
end
Base.propertynames(t::AbstractArray{<:Trace}, private=false) = fieldnames(eltype(t))
Base.:(==)(t1::Trace, t2::Trace) =
    all(x -> isequal(x[1], x[2]), (getfield.((t1, t2), f) for f in TRACE_FIELDS))

# Treat single objects as scalars in broadcasting
Base.broadcastable(t::Union{AbstractTrace,Event,Station,Position}) = Ref(t)

# Element type of trace
Base.eltype(::Trace{T,V,P}) where {T,V,P} = eltype(V)
