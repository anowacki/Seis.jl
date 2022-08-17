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

julia> Geodesy.ENU(evt::GeogEvent) = ENU(evt.lat, evt.lon, evt.elev)
```
"""
const GeogEvent{T} = Event{T, Geographic{T}}

const EVENT_FIELDS = fieldnames(Event)

function Base.getproperty(e::AbstractArray{<:Event}, f::Symbol)
    if f === :lon || f === :lat || f === :dep || f === :x || f === :y || f === :z
        getproperty.(e, f)
    elseif f === :time || f === :meta || f === :id
        getproperty.(e, f)
    else
        getfield(e, f)
    end
end

Base.setproperty!(e::AbstractArray{<:Event}, f::Symbol, val) =
    setproperty!.(e, f, val)
Base.propertynames(e::AbstractArray{<:Event}, private=false) = fieldnames(eltype(e))
Base.:(==)(e1::Event, e2::Event) =
    all(x -> isequal(x[1], x[2]), (getfield.((e1, e2), f) for f in EVENT_FIELDS))
