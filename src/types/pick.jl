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
