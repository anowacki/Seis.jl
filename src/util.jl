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
    ismissing(t.evt.time) && error("trace does not have origin time set")
    t.delta < 1e-3 && error("dates does not support sampling intervals < 1 ms")
    delta = Millisecond(round(Int, t.delta*1e3))
    t.evt.time:delta:(t.evt.time + (nsamples(t)-1)*delta)
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
