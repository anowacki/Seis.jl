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
Base.propertynames(s::AbstractArray{<:Station}, private=false) = propertynames(first(s))
Base.:(==)(s1::Station, s2::Station) =
    all(x -> isequal(x[1], x[2]), (getfield.((s1, s2), f) for f in STATION_FIELDS))
