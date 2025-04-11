# Types used by Seis

"Default floating point type used for traces if none specified as type parameters"
const DEFAULT_FLOAT = Float64

"""
    SeisDict

Wrapper around `Base.Dict` which allows one to get or set values using
the `{get|set}property[!]` syntax, i.e., `.`-access like `dict.key`.
`SeisDict`s also differ in that `missing` is returned instead of throwing a
`KeyError` when accessing a nonexistent key.  A key is removed if it is set to
`missing`.

# Examples

```
julia> dict = Seis.SeisDict(:a=>1)
Seis.SeisDict{Any,Any} with 1 entry:
  :a => 1

julia> dict.a
1

julia> dict.b
missing

julia> dict.a = missing
missing

julia> dict
Seis.SeisDict{Any,Any} with 0 entries
```

Access via `getindex` and `setindex!` (using `[]`s) is still possible:

```
julia> dict[:c] = 3
3

julia> dict[:c]
3

julia> dict[:c] = missing
missing

julia> dict[:d] # No key with this value
missing

julia> dict
Seis.SeisDict{Any,Any} with 0 entries
```

Arrays of `SeisDict`s also define `.`-access and setting via broadcasting, so one
may do:

```
julia> d1 = Seis.SeisDict(:a=>1)
Seis.SeisDict{Any,Any} with 1 entry:
  :a => 1

julia> d2 = deepcopy(d1); d2.a = 2;

julia> d = [d1, d2]
2-element Array{Seis.SeisDict{Any,Any},1}:
 Seis.SeisDict(:a => 1)
 Seis.SeisDict(:a => 2)

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
Base.setproperty!(sd::SeisDict{K,V}, key::Symbol, ::Missing) where {K,V} = delete!(sd, key)
Base.setproperty!(sd::SeisDict{K,Any}, key::Symbol, ::Missing) where K = delete!(sd, key)
Base.setproperty!(sd::SeisDict{K,V}, key::Symbol, val) where {K,V} = setindex!(sd, V(val), key)
Base.setproperty!(sd::SeisDict{K,Any}, key::Symbol, val) where K = setindex!(sd, val, key)
Base.propertynames(sd::SeisDict, private::Bool=false) = collect(keys(sd))

function Base.getproperty(sd::AbstractArray{<:SeisDict}, key::Symbol)
    if key in fieldnames(typeof(sd))
        getfield(sd, key)
    else
        getproperty.(sd, key)
    end
end
function Base.setproperty!(sd::AbstractArray{<:SeisDict}, key::Symbol, val)
    setproperty!.(sd, key, val)
end
# Method for concrete type needed in Julia v1.11+ now that Array is a Julia
# type with method.  This is type piracy!
function Base.setproperty!(sd::Array{<:SeisDict}, key::Symbol, val)
    setproperty!.(sd, key, val)
end

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

function Base.getproperty(p::AbstractArray{<:Position}, f::Symbol)
    if f in fieldnames(typeof(p))
        getfield(f, p)
    else
        getproperty.(p, f)
    end
end
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
- `elev`: Elevation above the reference level (e.g., ellipsoid) (m)
- `dep`: Depth below the reference level (km) (defined to be `-elev/1000`)

!!! note
    The fields `.elev` and `.dep` represent the same information, but in
    different ways.  Setting `dep` will automatically update `elev` and
    vice versa.  Both are provided for convenience.
"""
mutable struct Geographic{T<:AbstractFloat} <: Position{T}
    lon::Union{Missing,T}
    lat::Union{Missing,T}
    elev::Union{Missing,T}
end

"""
    Geographic(lon, lat, elev)
    Geographic{T}(lon, lat, elev)

Constructors for `Geographic` using positional arguments.
If passed with no type parameter `{T}` then by default `$DEFAULT_FLOAT` is used.
"""
Geographic(args...; kwargs...) = Geographic{DEFAULT_FLOAT}(args...; kwargs...)

"""
    Geographic(; lon, lat, elev, dep)
    Geographic{T}(; lon, lat, elev, dep)

Keyword argument constructors; the default for all values is `missing`.

!!! note
    If both `elev` and `dep` are passed, then an error is thrown if
    they do not equate to the same elevation above the reference level
    (i.e., if `dep != -elev/1000`).
"""
function Geographic{T}(; lon=missing, lat=missing, elev=missing, dep=missing) where T
    if dep !== missing && elev !== missing && dep != -elev/1000
        throw(ArgumentError(
            "cannot set both dep and elev unless they represent the same elevation"
        ))
    elseif dep !== missing
        elev = -1000*dep
    end
    Geographic{T}(lon, lat, elev)
end

Base.propertynames(::Geographic) = (:lon, :lat, :elev, :dep)

function Base.getproperty(g::Geographic, p::Symbol)
    if p === :dep
        -getfield(g, :elev)/1000
    else
        getfield(g, p)
    end
end

function Base.setproperty!(g::Geographic{T}, p::Symbol, v) where T
    if p === :dep
        setfield!(g, :elev, convert(Union{Missing,T}, -1000*v))
    else
        setfield!(g, p, convert(Union{Missing,T}, v))
    end
end

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

"""
    Cartesian(x, y, z)
    Cartesian{T}(x, y, z)

Constructors for `Cartesian` using positional arguments.
If passed with no type parameter `{T}` then by default `$DEFAULT_FLOAT` is used.
"""
Cartesian(args...; kwargs...) = Cartesian{DEFAULT_FLOAT}(args...; kwargs...)

"""
    Cartesian(; x, y, z)
    Cartesian{T}(; x, y, z)

Keyword argument constructors; the default for all values is `missing`.
"""
Cartesian{T}(; x=missing, y=missing, z=missing) where T =
    Cartesian{T}(x, y, z)

# Fix getproperty(::Array{<:Position}, ::Symbol) for Julia v1.11+
for T in (Cartesian, Geographic)
    @eval function Base.getproperty(a::Array{<:$T}, f::Symbol)
        if f in fieldnames(typeof(a))
            getfield(a, f)
        else
            getproperty.(a, f)
        end
    end
end

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
    Event{T,P}(P(lon, lat, -1000*dep), time, id, meta)

Event{T,P}(; x=missing, y=missing, z=missing, time=missing,
        id=missing, meta=SeisDict{Symbol,Any}()) where {T, P<:Cartesian} =
    Event{T,P}(P(x, y, z), time, id, meta)

Event{T}(; kwargs...) where T = Event{T ,Geographic{T}}(; kwargs...)
Event(; kwargs...) = Event{DEFAULT_FLOAT, Geographic{DEFAULT_FLOAT}}(; kwargs...)

function Base.propertynames(e::Event, private::Bool=false)
    pos_fields = filter(x->x!==:elev, propertynames(getfield(e, :pos)))
    if private
        (pos_fields..., fieldnames(Event)...)
    else
        (pos_fields..., filter(x->x!==:pos, fieldnames(Event))...)
    end
end

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

julia> Geodesy.ENU(evt::GeogEvent) = ENU(evt.lat, evt.lon, evt.elev)
```
"""
const GeogEvent{T} = Event{T, Geographic{T}}

const EVENT_FIELDS = fieldnames(Event)

function Base.getproperty(e::AbstractArray{<:Event}, f::Symbol)
    if f === :lon || f === :lat || f === :dep || f === :elev || f === :x || f === :y || f === :z
        getproperty.(e, f)
    elseif f === :time || f === :meta || f === :id
        getproperty.(e, f)
    else
        getfield(e, f)
    end
end

Base.setproperty!(e::AbstractArray{<:Event}, f::Symbol, val) =
    setproperty!.(e, f, val)
Base.setproperty!(e::Array{<:Event}, f::Symbol, val) = setproperty!.(e, f, val)
Base.propertynames(e::AbstractArray{<:Event}, private=false) = fieldnames(eltype(e))
Base.:(==)(e1::Event, e2::Event) =
    all(x -> isequal(x[1], x[2]), (getfield.((e1, e2), f) for f in EVENT_FIELDS))

"""
    Station

Struct containing information about a seismic station.  Fields `net`, `sta`, `loc`
and `cha` are the station, network, channel and location codes respectively, whilst
`lon` and `lat` are the location in °.  Set station elevation `elev`
in m relative to the reference level.  The azimuth `azi` and inclination `inc` of
the channel in ° are respectively measured from north to east, and downward from the
vertical.  (E.g., a "BHN" channel typically will have a `azi == 0` and `inc == 90`.)

`meta` is a `Dict` holding any extra
information about the station or channel, such as burial depth in m.

Missing information is allowed and stored as `missing`.
"""
mutable struct Station{T<:AbstractFloat, P<:Position{T}}
    net::Union{String,Missing}
    sta::Union{String,Missing}
    loc::Union{String,Missing}
    cha::Union{String,Missing}
    pos::P
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
- `azi`: component azimuth, in ° east of local north
- `inc`: component inclination, in ° down from upwards
- `meta`: `Dict` holding any extra information

By default, stations are geographic and `P <: Seis.Geographic`, in which
case the following may be passed:
- `lon`: longitude in °
- `lat`: latitude in °
- `elev`: elevation above reference level (usually the ellipsoid) in m

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
        cha=missing, lon=missing, lat=missing, elev=missing,
        azi=missing, inc=missing, meta=Dict()) where {T, P<:Geographic} =
    Station{T, P}(net, sta, loc, cha, P(lon, lat, elev), azi, inc, meta)

Station{T,P}(; net=missing, sta=missing, loc=missing,
        cha=missing, x=missing, y=missing, z=missing, elev=missing,
        azi=missing, inc=missing, meta=Dict()) where {T, P<:Cartesian} =
    Station{T,P}(net, sta, loc, cha, P(x, y, z), azi, inc, meta)

Station{T}(; kwargs...) where T = Station{T, Geographic{T}}(; kwargs...)
Station(; kwargs...) = Station{DEFAULT_FLOAT, Geographic{DEFAULT_FLOAT}}(; kwargs...)

function Base.propertynames(s::Station, private::Bool=false)
    pos_fields = filter(x->x!==:dep, propertynames(getfield(s, :pos)))
    if private
        (pos_fields..., fieldnames(Station)...)
    else
        (pos_fields..., filter(x->x!==:pos, fieldnames(Station))...)
    end
end

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
    if p === :lon || p === :lat || p === :dep || p === :elev || p === :x || p === :y || p === :z
        getproperty(getfield(e, :pos), p)
    else
        getfield(e, p)
    end
end

function Base.setproperty!(e::Union{Event{T},Station{T}}, p::Symbol, v) where T
    if p === :lon || p === :lat || p === :dep || p === :elev || p === :x || p === :y || p === :z
        setproperty!(getfield(e, :pos), p, convert(Union{Missing,T}, v))
    else
        setfield!(e, p, convert(fieldtype(typeof(e), p), v))
    end
end

function Base.getproperty(s::AbstractArray{<:Station}, f::Symbol)
    if f === :lon || f === :lat || f === :dep || f === :elev
        getproperty.(s, f)
    elseif f === :x || f === :y || f === :z
        getproperty.(s, f)
    elseif f === :net || f === :sta || f === :cha || f == :loc
        getproperty.(s, f)
    elseif f === :inc || f === :azi || f === :meta
        getproperty.(s, f)
    else
        getfield(s, f)
    end
end

Base.setproperty!(s::AbstractArray{<:Station}, f::Symbol, val) =
    setproperty!.(s, f, val)
Base.setproperty!(s::Array{<:Station}, f::Symbol, val) = setproperty!.(s, f, val)
Base.propertynames(s::AbstractArray{<:Station}, private=false) = propertynames(first(s))
Base.:(==)(s1::Station, s2::Station) =
    all(x -> isequal(x[1], x[2]), (getfield.((s1, s2), f) for f in STATION_FIELDS))

"""
    Pick{T}

A named arrival time which can be associated with a `Trace`.

If `p` is a `Pick`, then `p.time` gives the arrival time, and `p.name` gives its
description.

`Pick`s can be iterated to retrieve the time and name in that order.
Collections of `Pick`s (such as those stored in the `picks` field of `Trace`s)
can have their elements assigned by automatic conversion from real numbers,
tuples, or named tuples with fields `time` and `name`.

Arrays of `Pick`s, which are stored in `Traces`, can be accessed using
`getproperty` (`.`-access) and this returns arrays of times or names.

# Examples

## Getting and setting values

```
julia> t = Trace(0, 1, 1)
Seis.Trace{Float64,Array{Float64,1},Seis.Geographic{Float64}}:
            b: 0.0
        delta: 1.0
 Station{Float64,Seis.Geographic{Float64}}:
     sta.meta: Seis.SeisDict{Symbol,Any}()
 Event{Float64,Seis.Geographic{Float64}}:
     evt.meta: Seis.SeisDict{Symbol,Any}()
 Trace:
        picks: 0
         meta: 

julia> t.picks.P = 1.2
1.2

julia> t.picks.S = 2.1, "Sn"
(2.1, "Sn")

julia> t.picks.S.time
2.1

julia> t.picks
Seis.SeisDict{Union{Int64, Symbol},Seis.Pick{Float64}} with 2 entries:
  :P => Seis.Pick{Float64}(time=1.2, name=missing)
  :S => Seis.Pick{Float64}(time=2.1, name="Sn")
```

## Iteration

```
julia> time, name = Seis.Pick(1, "X")
Seis.Pick{Float64}(time=1.0, name="X")

julia> time, name
(1.0, "X")
```
"""
struct Pick{T<:AbstractFloat}
    time::T
    name::Union{Missing,String}
    Pick{T}(time, name=missing) where T = new{T}(time, name)
end

Pick(time, name) = Pick{DEFAULT_FLOAT}(time, name)
Pick(args) = Pick{DEFAULT_FLOAT}(args...)
Pick(; kwargs...) = Pick{DEFAULT_FLOAT}(; kwargs...)
Pick{T}(tup::Tuple{<:Any,<:Any}) where T = Pick{T}(tup...)
Pick{T}(; time, name=missing) where T = Pick{T}(time, name)
Base.convert(::Type{Pick}, time::Real) = Pick(time)
Base.convert(::Type{Pick}, x::Tuple{<:Any,<:Any}) = Pick(x...)
Base.convert(::Type{Pick}, p::Pick) = Pick(p.time, p.name)
Base.convert(::Type{Pick}, x::NamedTuple{(:time,:name)}) = Pick(x...)
Base.convert(::Type{Pick{T}}, time::Real) where T = Pick{T}(time)
Base.convert(::Type{Pick{T}}, x::Tuple{<:Any,<:Any}) where T = Pick{T}(x...)
Base.convert(::Type{Pick{T}}, p::Pick) where T = Pick{T}(p.time, p.name)
Base.convert(::Type{Pick{T}}, x::NamedTuple{(:time,:name)}) where T = Pick{T}(x...)
Base.first(p::Pick) = p.time
Base.last(p::Pick) = p.name

# Picks iterate time, name
Base.iterate(p::Pick, i=1) = i == 1 ? (p.time, 2) : (i == 2 ? (p.name, 3) : nothing)

Base.getindex(p::Pick, i::Integer) = getfield(p, i)

# Vectors of picks can have their values retrieved as a vector
function Base.getproperty(p::AbstractArray{<:Pick}, field::Symbol)
    if field === :time || field === :name
        getproperty.(p, field)
    else
        getfield(p, field)
    end
end

"""
    AbstractData

Abstract supertype of all single-channel datatypes in Seis.jl.  An
`AbstractData` represents some data acquired at a point where it
makes sense to associate that recording with a channel, and where the
recording is limited to some period of time.  However, the recording
may be in the time or frequency domain, and need not be stationary or
evenly-sampled.
"""
abstract type AbstractData end

"""
    AbstractTrace <: AbstractData

Abstract type from which you should subtype if creating new types of traces.
`AbstractTrace`s are time-domain recordings (or synthetics) at a single
channel.

# Interface

!!! note
    The formal interface for `AbstractTrace`s is still a work in progress and may change with
    a minor version increment.

The following methods should be defined for all `AbstractTrace`s `t`:

- `trace(t)`: Return the data for the trace.
- `times(t)`: Return the time at each sample of `t`.
- `starttime(t)`: The time of the first sample.
- `nsamples(t)`: The number of samples in `t`.
- `Base.eltype(t)`: The element type of the data samples.
- `t.evt`: Return the `Event` associated with this trace.
- `t.sta`: Return the `Station` at which this trace was recorded.
- `t.meta`: Return a `SeisDict{Symbol,Any}` into which metadata may be placed.
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

---

    Trace(b, delta, n::Integer) -> trace::Trace{$DEFAULT_FLOAT,Vector{$DEFAULT_FLOAT},Seis.Geographic{$DEFAULT_FLOAT}}

Create a new `Trace` with uninitialised data of length `n` samples.

---

    Trace(; kwargs...) -> trace

Construct a new `Trace` with fields set according the to the following
keyword arguments `kwargs`:

- `b = 0`: Start time of trace relative to `evt.time`
- `delta = 1`: Sampling interval in s
- `data = $(DEFAULT_FLOAT)[]`: `AbstractVector` of trace data.  Takes
  precedence over `n` if `n` is also supplied.
- `n`: Create a vector with arbitrary contents of length `n` samples.
  Ignored if `data` is supplied.
- `evt = Event{$(DEFAULT_FLOAT)}()`: [`Event`](@ref) object
- `sta = Station{$(DEFAULT_FLOAT)}()`: [`Station`](@ref) object
- `meta = Seis.SeisDict{Symbol,Any}()`: Metadata dictionary
- `picks = Seis.SeisDict{Union{Int,Symbol},Seis.Pick{T}}()`: Picks dictionary

---

    Trace{T,V,P}(args...; kwargs...) -> trace::Trace{T,V,P}
    Trace{T,V}(args...; kwargs...) -> trace::Trace{T,V,Seis.Geographic{T}}
    Trace{T}(args...; kwargs...) -> trace::Trace{T,Vector{T},Seis.Geographic{T}}

Construct traces with non-default number and data types.  In the third
form, the data type defaults to `Vector{$DEFAULT_FLOAT}`, whilst in
both the second and third forms, `P` defaults to `Seis.Geographic{T}`.

# Examples
The default, empty trace
```
julia> Trace()
```

A trace with defined station code
"""
Trace{T,V,P}(b, delta, t::AbstractVector) where {T,V,P} =
    Trace{T,V,P}(b, delta, t, Event{T,P}(), Station{T,P}(), Dict(), Dict())
Trace{T,V,P}(b, delta, n::Integer) where {T,V,P} = Trace{T,V,P}(b, delta, Vector{T}(undef, n))
Trace{T,V}(args...) where {T,V} = Trace{T,V,Geographic{T}}(args...)
Trace{T}(args...) where T = Trace{T,Vector{T},Geographic{T}}(args...)
Trace(b, delta, t_or_n) = Trace{DEFAULT_FLOAT,Vector{DEFAULT_FLOAT},Geographic{DEFAULT_FLOAT}}(b, delta, t_or_n)

function Trace{T,V,P}(;
    b=zero(T),
    delta=one(T),
    n=0,
    data=V(undef, n),
    sta=Station{T,P}(),
    evt=Event{T,P}(),
    meta=SeisDict{Symbol,Any}(),
    picks=SeisDict{Union{Int,Symbol},Pick{T}}()
) where {T,V,P}
    if length(data) != n && n != 0
        @warn("ignoring keyword argument `n` as `data` was passed in")
    end
    Trace{T,V,P}(b, delta, data, evt, sta, picks, meta)
end
Trace{T,V}(; kwargs...) where {T,V} = Trace{T,V,Geographic{T}}(; kwargs...)
Trace{T}(; kwargs...) where T = Trace{T,Vector{T},Geographic{T}}(; kwargs...)
Trace(; kwargs...) = Trace{DEFAULT_FLOAT,Vector{DEFAULT_FLOAT},Geographic{DEFAULT_FLOAT}}(; kwargs...)

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
    CartTrace{T}(args...; kwargs...)

Construct a `Trace` with Cartesian coordinates for the `Event` and `Station`
with non-default number and data types.  See [`Trace{T,V}`](@ref)
for details.
"""
CartTrace{T}(args...; kwargs...) where T = Trace{T, Vector{T}, Cartesian{T}}(args...; kwargs...)
CartTrace(args...; kwargs...) = Trace{DEFAULT_FLOAT, Vector{DEFAULT_FLOAT},
          Cartesian{DEFAULT_FLOAT}}(args...; kwargs...)

const TRACE_FIELDS = fieldnames(Trace)

# Get an array of values from an array of traces
function Base.getproperty(t::AbstractArray{<:Trace}, f::Symbol)
    if f === :b || f === :delta || f === :meta || f === :evt || f === :sta || f === :picks
        getfield.(t, f)
    else
        getfield(t, f)
    end
end

for A in (AbstractArray{<:Trace}, Array{<:Trace})
    @eval begin
        # Set an array of traces with a single value
        function Base.setproperty!(t::$A, f::Symbol, val)
            if f === :b || f === :delta || f === :meta || f === :evt || f == :sta || f === :picks
                for tt in t
                    setproperty!(tt, f, val)
                end
            else
                # Fallback in case of desired dot-access to fields of an array type
                setfield!(t, f, val)
            end
        end

        # Set an array of traces with an array of values
        function Base.setproperty!(t::$A, f::Symbol, val::AbstractArray)
            length(t) == length(val) || throw(DimensionMismatch())
            for (tt, vv) in zip(t, val)
                setproperty!(tt, f, vv)
            end
        end
    end
end


Base.propertynames(t::AbstractArray{<:Trace}, private=false) = fieldnames(eltype(t))

Base.:(==)(t1::Trace, t2::Trace) =
    all(x -> isequal(x[1], x[2]), (getfield.((t1, t2), f) for f in TRACE_FIELDS))

# Define hash so that Base.isequal (and so Base.unique) work correctly
# for our types.  SeisDict does the right thing already.
for T in (Cartesian, Geographic, Event, Station, Trace)
    Tsym = nameof(T)
    fields = fieldnames(T)
    @eval function Base.hash(x::$Tsym, h::UInt)
        $([:(h = hash(x.$f, h)) for f in fields]...)
        h
    end
end

# Treat single objects as scalars in broadcasting
Base.broadcastable(t::Union{Trace,Event,Station,Position,Pick}) = Ref(t)

# Element type of trace
Base.eltype(::Trace{T,V,P}) where {T,V,P} = eltype(V)
