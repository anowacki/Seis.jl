"""
    angle_difference(α, β, degrees::Bool=true) -> Δ

Return the angular difference `Δ` between two angles `α` and `β`, in the direction
α → β (i.e., Δ = (β - α)).  This means that Δ is positive if the direction β is
clockwise of α, and negative otherwise.

Angles are assumed to be degreees unless `degrees` is `false`.
"""
function angle_difference(a::Real, b::Real, degrees::Bool=true)
    whole_circ = degrees ? 360.0 : 2π
    half_circ = whole_circ/2
    mod(b - a + half_circ, whole_circ) - half_circ
end

"""
    dates(t) -> date_range

Return a `date_range` which contains the dates for each sample of `t`, so long
as `t.evt.time` is defined.  If not, an error is thrown.

N.B.  This function assumes that the sampling interval `t.delta` is representable
as an integer number of milliseconds, and rounds it accordingly.  `Dates.DateTime`s
have precision of 1 ms.  An error is thrown if `t.delta < 1e-3` s.
"""
function dates(t)
    b, delta = _check_date_b_delta(t)
    (t.evt.time + b):delta:(t.evt.time + b + (nsamples(t)-1)*delta)
end

"""
    startdate(t) -> date

Return the `date` of the first sample of the trace `t`.

N.B.  This function assumes that the sampling interval `t.delta` is representable
as an integer number of milliseconds, and rounds it accordingly.  `Dates.DateTime`s
have precision of 1 ms.  An error is thrown if `t.delta < 1e-3` s.
"""
startdate(t::AbstractTrace) = ((b, delta) = _check_date_b_delta(t); t.evt.time + b)

"""
    enddate(t) -> date

Return the `date` of the last sample of the trace `t`.

N.B.  This function assumes that the sampling interval `t.delta` is representable
as an integer number of milliseconds, and rounds it accordingly.  `Dates.DateTime`s
have precision of 1 ms.  An error is thrown if `t.delta < 1e-3` s.
"""
enddate(t::AbstractTrace) = ((b, delta) = _check_date_b_delta(t); t.evt.time + b + (nsamples(t)-1)*delta)

"""
    _check_date_b_delta(t) -> b::Millisecond, delta::Millisecond

Throw an error if a trace either has no origin time set, or has a sampling interval
less than 1 ms, and return the trace start time `b` and sampling interval `delta`
as `Dates.Millisecond`s.
"""
function _check_date_b_delta(t::AbstractTrace)
    ismissing(t.evt.time) && error("trace does not have origin time set")
    t.delta < 1e-3 && error("date calculations do not support sampling intervals < 1 ms")
    b = Millisecond(round(Int, starttime(t)*1000))
    delta_in_ms = round(Int, 1000*t.delta)
    delta_in_ms ≈ 1000t.delta || error("date calculations do not support sampling intervals " *
                                       "which are not whole numbers of milliseconds")
    delta = Millisecond(delta_in_ms)
    b, delta
end

"""
    starttime(t) -> time

Return the start `time` of trace `t` in seconds.
"""
starttime(t::AbstractTrace) = t.b

"""
    endtime(t) -> time

Return the end `time` of trace `t` in seconds.
"""
endtime(t::AbstractTrace) = t.b + (nsamples(t) - 1)*t.delta

"""
    origin_time!(t, time::DateTime; picks=true) -> t

Set the origin time of the trace `t` and shift the start time of the trace
(stored in its `.b` field) so that the absolute time of all samples remains
the same.

`origin_time!` will also shift all pick times so that they remain at
the same absolute time.  Set `picks=false` to leave picks at the same
time relative to the trace start time.

If `t.evt.time` is `missing` (i.e., unset), then it is simply set to
`time` and no times are shifted.
"""
function origin_time!(t::AbstractTrace, time::DateTime; picks=true)
    if t.evt.time === missing
        t.evt.time = time
        return t
    end
    # Shift in s of the trace start time from old to new origin time.
    # Δb is positive if the new time is *later*.
    Δb = Dates.value(time - t.evt.time)*1e-3
    t.b -= Δb
    t.evt.time = time
    if picks
        for (key, (time, name)) in t.picks
            t.picks[key] = (time=time-Δb, name=name)
        end
    end
    t
end

"""
    origin_time(t, time::DateTime; picks=true) -> t′

Return a copy to `t` where the event origin time is shifted to `time`.

See the in-place version [`origin_time!`](@ref) for more details.
"""
origin_time(t::AbstractTrace, time::DateTime; kwargs...) =
    origin_time!(deepcopy(t), time; kwargs...)

"""
    is_horizontal(s::Station{T}; tol=eps(T)) where T
    is_horizontal(t::AbstractTrace; tol=eps(eltype(trace(t)))) -> ::Bool

Return `true` if the trace `t` is horizontal (i.e., its inclination
is 90° from the vertical), and `false` otherwise.
"""
is_horizontal(s::Station{T}; tol=eps(T))  where T = isapprox(s.inc, 90, atol=tol)
is_horizontal(t::Trace{T}; tol=eps(T)) where T = is_horizontal(t.sta, tol=tol)

"""
    is_vertical(s::Station{T}; tol=eps(T)) where T
    is_vertical(t::AbstractTrace; tol=eps(eltype(trace(t)))) -> ::Bool

Return `true` if the trace `t` is vertical (i.e., its inclination
is 0°), and `false` otherwise.
"""
is_vertical(s::Station{T}; tol=eps(T))  where T = isapprox(s.inc, 0, atol=tol)
is_vertical(t::Trace{T}; tol=eps(T)) where T = is_vertical(t.sta, tol=tol)

"""
    linear_regression(x, y)

Perform simple linear regression using Ordinary Least Squares. Returns `a` and `b` such
that `a + b*x` is the closest straight line to the given points `(x, y)`, i.e., such that
the squared error between `y` and `a + b*x` is minimized.
"""
function linear_regression(x::AbstractVector, y::AbstractVector)
    size(x) == size(y) || throw(DimensionMismatch("x and y must be the same size"))
    mx, my = mean(x), mean(y)
    b = covm(x, mx, y, my)/varm(x, mx)
    a = my - b*mx
    return a, b
end

"""
    nearest_sample(t::AbstractTrace, time; inside=true) -> i

Return the index `i` of the nearest sample of the trace `t` to `time` seconds.

If `inside` is `true` (the default), return `nothing` when `time` lies outside
the trace.  Set `inside` to `false` to instead return the first or last index
when `time` is outside the trace.
"""
function nearest_sample(t::AbstractTrace, time; inside=true)::Union{Int,Nothing}
    if inside
        t.b <= time <= endtime(t) || return nothing
    end
    time <= t.b && return 1
    time >= endtime(t) && return nsamples(t)
    round(Int, (time - t.b)/t.delta + 1)
end

"""
    nearest_sample(t::AbstractTrace, datetime::DateTime; inside=true)

Form of `nearest_sample` where `datetime` is given as absolute time.

An error is thrown if no origin time is specified for `t.evt.time`.
"""
function nearest_sample(t::AbstractTrace, datetime::DateTime; inside=true)
    ismissing(t.evt.time) && error("trace does not have origin time set")
    d = dates(t)
    if inside
        t.evt.time <= datetime <= d[end] || return nothing
    end
    datetime <= t.evt.time && return 1
    datetime >= d[end] && return nsamples(t)
    argmin(abs.(d .- datetime))
end

"""
    nsamples(t) -> n

Return the number of samples `n` in a trace `t`.
"""
nsamples(t::AbstractTrace)::Int = length(t.t)

"""
    nsamples(t, b, e) -> n

Return the number of samples `n` in a trace `t` between times `b`
and `e` seconds.

This function only counts samples that are strictly on or later than
the `b` time, and before or on the `e` time.

# Example
```
julia> t = Trace(0, 1, 25); # 25 samples from 0s to 24 s

julia> nsamples(t, 3, 4.1)
2
```

See also: [`nearest_sample`](@ref).
"""
function nsamples(t::AbstractTrace, b, e)
    n = nsamples(t)
    tb = starttime(t)
    te = endtime(t)
    e < tb && return 0
    b > te && return 0
    delta = t.delta
    ib = clamp(ceil(Int, (b - tb)/delta) + 1, 1, n)
    ie = clamp(floor(Int, (e - tb)/delta) + 1, 1, n)
    max(0, ie - ib + 1)
end

"""
    nsamples(t, start::DateTime, stop::DateTime) -> n

Return the number of samples `n` in a trace `t` between dates
`start` and `stop`.

# Example
```
julia> using Dates

julia> t = Trace(10, 1, 20); # 20 samples from 10 to 30 s

julia> t.evt.time = DateTime(3000)
3000-01-01T00:00:00

julia> nsamples(t, DateTime(3000), DateTime(3000) + Second(9))
0

julia> nsamples(t, DateTime(3000) + Second(20), DateTime(3000) + Second(22))
3
```
"""
function nsamples(t::AbstractTrace, start::DateTime, stop::DateTime)
    ismissing(t.evt.time) && throw(ArgumentError("trace does not have origin time set"))
    bdate = startdate(t)
    edate = enddate(t)
    (start > edate || stop < bdate) && return 0
    tb = starttime(t)
    b = Dates.value(start - bdate)/1000 + tb
    e = Dates.value(stop - bdate)/1000 + tb
    nsamples(t, b, e)
end

"""
    times(t) -> range

Return the set of times `range` at which each sample of the trace `t` is defined.
"""
times(t::AbstractTrace) = t.b .+ (0:(nsamples(t) - 1)).*t.delta

"""
    trace(t) -> y

Return an array containing the values of the `Trace` `t`.
"""
trace(t::AbstractTrace) = t.t

"""
    traces_are_orthogonal(t1::Trace, t2::Trace; tol=eps()) -> ::Bool

Return `true` if the two traces `t1` and `t2` have component azimuths
90° apart.  Set the tolerance of the comparison with `tol`.
"""
traces_are_orthogonal(t1::AbstractTrace, t2::AbstractTrace;
                      tol=max(√eps(eltype(trace(t1))), √eps(eltype(trace(t2))))) =
    isapprox(abs(angle_difference(t1.sta.azi, t2.sta.azi)), 90, atol=tol)

"""
    @chain function f(t::Trace, args...; kwargs...) ... end
    @chain f(t::Trace, args...; kwargs...) = ...

Allow function `f` to be used in chaining where the first argument is a `Trace`.
This allows, for instance:

```
julia> using Seis

julia> Seis.@chain b_offset(t::Trace, offset) = t.b + offset
b_offset (generic function with 1 method)

julia> Trace(0, 1, [1]) |> b_offset(2)
2.0
```

The chaining functions that Seis offers are implemented using this macro.

### Known limitations

Currently the @chain macro doesn't cope with docstrings, which must be added
separately like:

```
@chain bandpass(t::Trace, low, high) = ...

\"Documentation for bandpass\"
bandpass
```
"""
macro chain(ex)
    chain_func = if @capture(ex, f_(t_::Trace, args__; kwargs__) = body_) ||
            @capture(ex, function f_(t_::Trace, args__; kwargs__) body_ end)
        :($f($(args...); $(kwargs...)) = x -> $f(x, $(args...); $(kwargs...)))
    elseif @capture(ex, f_(t_::Trace, args__) = body_) ||
            @capture(ex, function f_(t_::Trace, args__) body_ end)
        :($f($(args...)) = $t -> $f($t, $(args...)))
    else
        error("@chain must be used to annotate functions with a Trace " *
              "as the first argument")
    end
    return quote
        $(ex)
        $(chain_func)
    end |> esc
end
