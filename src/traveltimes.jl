# Functions relating to travel times and picks

"""
    add_pick!(t, time [, name=missing]) -> ::Seis.Pick

Add an arrival time pick to the Trace `t`, ensuring existing picks are not
overwritten, and return the `Seis.Pick` object added to the trace

If `name` is not `missing`, then the key of this pick will be `Symbol(name)`,
unless another pick with the same key already exists.  In that case, the name
will be appended with a number which increases until an available key is found.

If `name` is missing, then the pick is added to a numbered set of picks.

(Direct manipulation of picks is easy: just do `t.picks.PKP = (1001, "PKP")` to set
a picks with name "PKP", time 1001 s and key `:PKP`.)

# Example

```
julia> t = Trace(0, 1, 2);

julia> add_pick!.(t, [1, 2], ["A", missing]);

julia> t.picks
Seis.SeisDict{Union{Int64, Symbol},Seis.Pick{Float64}} with 2 entries:
  :A => Seis.Pick{Float64}(time=1.0, name="A")
  1  => Seis.Pick{Float64}(time=2.0, name=missing)

julia> add_pick!(t, 4)
Seis.Pick{Float64}(time=4.0, name=missing)

julia> t.picks
Seis.SeisDict{Union{Int64, Symbol},Seis.Pick{Float64}} with 3 entries:
  :A => Seis.Pick{Float64}(time=1.0, name="A")
  2  => Seis.Pick{Float64}(time=4.0, name=missing)
  1  => Seis.Pick{Float64}(time=2.0, name=missing)

julia> t.picks.A
Seis.Pick{Float64}(time=1.0, name="A")

julia> t.picks[1]
Seis.Pick{Float64}(time=2.0, name=missing)
```

See also: [`Seis.Pick`](@ref).
"""
function add_pick!(t::AbstractTrace, time, name=missing)
    key = if !ismissing(name)
        base_key = name
        key = Symbol(base_key)
        i = 0
        while key in keys(t.picks)
            i += 1
            key = Symbol(base_key * "_" * string(i))
        end
        key
    else
        i = 1
        while i in keys(t.picks)
            i += 1
        end
        i
    end
    t.picks[key] = time, name
    t.picks[key]
end

"""
    add_pick!(t, date::Dates.AbstractDateTime[, name=missing]) -> ::Seis.Pick

Add a time pick based on absolute time, given as a `DateTime` or another
`AbstractDateTime`.

!!! note
    The pick is converted to relative time, so does not remain independent of
    any changes to `t.b` or `t.evt.time` if you update the trace manually.
    Use [`origin_time!`](@ref) to change the origin time whilst preserving
    the absolute time of picks.

# Example
```
julia> t = sample_data();

julia> using Dates

julia> add_pick!(t, DateTime("1981-03-29T10:39:10"), "Coffee time")
Seis.Pick{Float32}(time=56.0, name="Coffee time")

julia> picks(t)
3-element Vector{Seis.Pick{Float32}}:
 Seis.Pick{Float32}(time=53.670002, name=missing)
 Seis.Pick{Float32}(time=56.0, name="Coffee time")
 Seis.Pick{Float32}(time=60.980003, name=missing)
```
"""
function add_pick!(t::AbstractTrace, date::Dates.AbstractDateTime, name=missing)
    if origin_time(t) === missing
        throw(ArgumentError("origin time must be set for trace to add a pick by date"))
    end

    time_ms::Dates.Millisecond = date - origin_time(t)
    time = Dates.value(time_ms)/1000
    add_pick!(t, time, name)
end

"""
    add_pick!(t, p::Pick, name=p.name) -> p

Add a travel time pick to the `Trace` `t` from a `Seis.Pick`.  By default,
the pick name is used.

# Example
```
julia> t1 = Trace(0, 1, 20); t2 = sample_data();

julia> add_pick!(t1, t2.picks.A, "A")
Seis.Pick{Float64}(time=53.67000198364258, name="A")

julia> t1.picks
Seis.SeisDict{Union{Int64, Symbol},Seis.Pick{Float64}} with 1 entry:
  :A => Seis.Pick{Float64}(time=53.67000198364258, name="A")
```
"""
add_pick!(t::AbstractTrace, p::Pick, name=p.name) = add_pick!(t, p.time, name)

# Stub to allow methods to be added using a travel time package (like SeisTau)
"""
    add_picks!

Add picks to traces based on seismic phases' predicted arrival time.

Seis does not itself implement seismic phase travel time computation.
See [SeisTau](https://github.com/anowacki/SeisTau.jl) for one implementation
which adds methods to `add_pick!`.
"""
function add_picks! end

"""
    clear_picks!(t)

Remove all picks associated with the `Trace` `t`.

# Example
```
julia> t = sample_data();

julia> picks(t)
2-element Array{Seis.Pick{Float32},1}:
 Seis.Pick{Float32}(time=53.670002, name=missing)
 Seis.Pick{Float32}(time=60.980003, name=missing)

julia> clear_picks!(t);

julia> picks(t)
0-element Array{Seis.Pick{Float32},1}
```
"""
clear_picks!(t::AbstractTrace) = empty!(t.picks)

"""
    picks(t; sort=:time) -> p::Vector{Tuple{<:AbstractString,<:AbstractFloat}}

Return a vector `p` of `Seis.Pick`s, which contain pairs of pick times and names
associated with the `Trace` `t`.

`sort` can be one of:
- `:time` (the default): Picks are returned in order of increasing time
- `:name`: Picks are sorted alphanumerically by name, with unnamed picks first

The returned vector can be iterated like:

```
julia> t = Trace(0, 1, rand(10));

julia> add_pick!.(t, (1,2), ("P","S"));

julia> for (time, name) in picks(t) @show time, name end
(time, name) = (1.0, "P")
(time, name) = (2.0, "S")
```
"""
function picks(t::AbstractTrace; sort=:time)
    ps = collect(values(t.picks))
    _sortpicks(ps, sort)
end

"""
    picks(t, name::AbstractString; sort=:time) -> p
    picks(t, pattern::Regex; sort=:time) -> p

Return a vector `p` of pairs of pick names and times associated with the `Trace` `t`
which either are exactly `name` or match the regular expression `pattern`.

By default, picks are returned in order of increasing time.  Use `sort=:name`
to sort alphanumerically by name (where unnamed picks appear first).

# Example
```
julia> t = Trace(0, 1, 2); t.picks.P = (1, "Pn"); t.picks.S = (1.8, "S");

julia> t.picks
Seis.SeisDict{Union{Int64, Symbol},Seis.Pick{Float64}} with 2 entries:
  :P => Seis.Pick{Float64}(time=1.0, name="Pn")
  :S => Seis.Pick{Float64}(time=1.8, name="S")

julia> picks(t, "S")
1-element Array{Seis.Pick{Float64},1}:
 Seis.Pick{Float64}(time=1.8, name="S")

julia> picks(t, r"^P")
1-element Array{Seis.Pick{Float64},1}:
 Seis.Pick{Float64}(time=1.0, name="Pn")
```
"""
function picks(t::AbstractTrace, name_or_match; sort=:time)
    ps = _picks(t, name_or_match)
    _sortpicks(ps, sort)::Vector{Pick{eltype(t)}}
end

_picks(t::AbstractTrace, name::AbstractString) = filter(x->coalesce(x.name, "")==name, picks(t))
_picks(t::AbstractTrace, match::Regex) = filter(x->occursin(match, coalesce(x.name, "")), picks(t))
function _picks(t::AbstractTrace, key::Symbol)
    p = t.picks[key]
    p === missing ? Pick{eltype(t)}[] : [p]
end

_sortpicks(ps, ::Nothing) = ps
function _sortpicks(ps, sort)
    inds = if sort == :time
        sortperm(getfield.(ps, :time))
    elseif sort == :name
        sortperm(string.(replace(getfield.(ps, :name), missing=>"")))
    else
        throw(ArgumentError("`sort` can be only `:time` or `:name`"))
    end
    ps[inds]
end
