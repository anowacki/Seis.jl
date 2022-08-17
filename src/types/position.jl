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
