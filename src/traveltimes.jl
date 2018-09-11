# Functions relating to travel times and picks

"""
    add_pick!(t, time [, name=missing]) -> (time, name)

Add an arrival time pick to the Trace `t`.
"""
add_pick!(t::AbstractTrace, time, name=missing) = push!(t.picks, (time=time, name=name))

"""
    add_pick!(t, p::TauPy.Phase) -> (time, name)

Add a travel time pick to the `Trace` `t` from a `TauPy.Phase` arrival.
"""
add_pick!(t::AbstractTrace, p::TauPy.Phase) = add_pick!(t, p.time, p.name)

"""
    add_picks!(t, phase; model="iasp91", exact=false)

Add travel time picks to the trace `t` for the 1D Earth `model`, for arrivals
with a name matching `phase`.

If `exact` is `true`, only phases which are an exact match for `phase` will
be added.

Available models are: $(TauPy.available_models()).
"""
function add_picks!(t::AbstractTrace, phase::AbstractString="all"; model="iasp91", exact=false)
    _check_headers_taup(t)
    arrivals = travel_time(t, phase; model=model)
    for arr in arrivals
        exact && arr.name != phase && continue
        add_pick!(t, arr.time, arr.name)
    end
    picks
end

"""
    clear_picks!(t)

Remove all picks associated with the `Trace` `t`.
"""
clear_picks!(t::AbstractTrace) = t.picks = []

"""
    picks(t) -> p::Vector{Tuple{<:AbstractString,<:AbstractFloat}}

Return a vector `p` of pairs of pick times and names associated with the `Trace` `t`.
This can be iterated like:

```
julia> t = Trace(0, 1, rand(10));

julia> add_pick!.(t, (1,2), ("P","S"));

julia> for (time, name) in picks(t) @show time, name end
(time, name) = (1.0, "P")
(time, name) = (2.0, "S")

```
"""
picks(t::AbstractTrace) = t.picks

"""
    picks(t, name::AbstractString) -> p
    picks(t, pattern::Regex) -> p

Return a vector `p` of pairs of pick names and times associated with the `Trace` `t`
which either are exactly `name` or match the regular expression `pattern`.
"""
picks(t::AbstractTrace, name::AbstractString) = filter(x->coalesce(x[2], "")==name, picks(t))
picks(t::AbstractTrace, match::Regex) = filter(x->occursin(match, coalesce(x[2], "")), picks(t))

"""
    travel_time(t, phase="all"; model="iasp91") -> phases::Vector{TauPy.Phase}

Return travel times for the geometry specified in the `Trace` `t` for the 1D
Earth `model`, for arrivals with a name matching `phase`.

Available models are: $(TauPy.available_models()).
"""
function travel_time(t::AbstractTrace, phase::AbstractString; model="iasp91", kwargs...)
    _check_headers_taup(t)
    TauPy.travel_time(t.evt.dep, distance_deg(t), phase; model=model, kwargs...)
end

"""
    travel_time(t, model="iasp91") -> ::Vector{Vector{TauPy.Phase}}

Return the list of `TauPy.Phase`s associated with each travel time pick in `t`.
"""
travel_time(t::AbstractTrace; model::AbstractString="iasp91") =
    [travel_time(t, name; model=model) for (_, name) in picks(t) if name[1] in ("p", "s", "P", "S")]

"Throw an error if a Trace doesn't contain the right headers to call TauPy.travel_time."
_check_headers_taup(t::AbstractTrace) = any(ismissing,
    (t.evt.lon, t.evt.lat, t.evt.dep, t.sta.lon, t.sta.lat)) &&
    throw(ArgumentError("Insufficient information in trace to compute travel times." *
                        "  (Need: evt.lon, evt.lat, evt.dep, sta.lon and sta.lat)")) ||
    nothing
