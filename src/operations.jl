# Trace operations

"""
    cut!(t, start, end; warn=true) -> t

Cut a `Trace` `t` in place between `start` and `end`.  An error is thrown if either
`start` or `end` are `missing`.

By default, a warning is shown if cut times lie outside the trace; set `warn` to `false`
to turn this off.
"""
function cut!(t::Trace, b, e; warn=true)
    (b === missing || e === missing) && throw(ArgumentError("Start or end cut time is `missing`"))
    e < b && throw(ArgumentError("End cut time ($e) is before start cut ($b)"))
    if b < t.b
        warn && @warn("Beginning cut time $b is before start of trace.  Setting to $(t.b).")
        b = t.b
    end
    b > endtime(t) &&
        throw(ArgumentError("Beginning cut time $b is later than end of trace ($(endtime(t)))."))
    if e > endtime(t)
        warn && @warn("End cut time $e is after end of trace.  Setting to $(endtime(t)).")
        e = endtime(t)
    end
    e < t.b && throw(ArgumentError("End cut time $e is earlier than start of trace (t.b)."))
    ib = round(Int, (b - t.b)/t.delta) + 1
    ie = nsamples(t) - round(Int, (endtime(t) - e)/t.delta)
    t.t = t.t[ib:ie]
    t.b += (ib - 1)*t.delta
    t
end

"""
    cut!(t, start_date, end_date; kwargs...) -> t

Cut a `Trace` `t` in place between dates `start_date` and `end_date`.
"""
cut!(t::Trace, b::DateTime, e::DateTime; kwargs...) =
    cut(t, Dates.value(Dates.Nanosecond(b - t.evt.time))/1e9,
        Dates.value(Dates.Nanosecond(e - t.evt.time))/1e9; kwargs...)

"""
    cut!(t, pick1, offset1, pick2, offset; kwargs...) -> t
    cut!(t, pick, offset1, offset2; kwargs...) ->

Cut a trace `t` in place between `offset1` s after the first pick with name `pick1`
and `offset2` s after `pick2`.

In the second form, both offsets are relative to `pick`
"""
cut!(t::Trace, pick1, offset1, pick2, offset2; kwargs...) =
    cut!(t, first(picks(t, pick1)).time + offset1, first(picks(t, pick2)).time + offset2; kwargs...)
cut!(t::Trace, pick1, offset1, offset2; kwargs...) =
    cut!(t, pick1, offset1, pick1, offset2; kwargs...)

"""
    cut(t, start, end; kwargs...) -> t′
    cut(t, start_date, end_date; kwargs...) -> t′
    cut(t, pick1, offset1, pick2, offset2; kwargs...) -> t′
    cut(t, pick, offset1, offset2; kwargs...) -> t′

Return a copy of the trace `t` cut between `start` and `end` s relative to the event
origin.  You may also specify a `start_date` and `end_date`, or choose times
`offset1` and `offset2` s relative to `pick1` and `pick2` respectively.  Both offset
times may also be specified relative to one `pick`.
"""
cut(t::Trace, args...; kwargs...) = cut!(deepcopy(t), args...; kwargs...)

# TODO: Enable antialiasing
"""
    decimate!(t, n; antialias=true)
    
Decimate the trace `t` by removing all except every `n` points.  The sampling interval
is increased `n` times.

If `antialias` is `false`, then no antialiasing filtering is applied during decimation.
This means the decimated trace may contain spurious signals.
"""
function decimate!(t, n::Integer; antialias=true)
    antialias && error("Antialias filtering not implemented; use `antialias=false`")
    1 <= n || throw(ArgumentError("n must be grater than 0 (supplied $n)"))
    n == 1 && return t
    t.t = t.t[1:n:end]
    t.delta *= n
    t
end

"""
    decimate(t, n; antialias=true)

Decimate the trace `t` by removing all except every `n` points.  The sampling interval
is increased `n` times.

If `antialias` is `false`, then no antialiasing filtering is applied during decimation.
This means the decimated trace may contain spurious signals.
"""
decimate(t::Trace, n; antialias=true) = decimate!(deepcopy(t), n; antialias=antialias)

"""
    normalise!(t::Trace, val=1) -> t

Normalise the trace `t` so that its maximum absolute amplitude is `val`,
and return the trace.
"""
function normalise!(t::AbstractTrace, val=1)
    maxval = maximum(abs, trace(t))
    t.t .*= val/maxval
    t
end
normalise(t::AbstractTrace, args...; kwargs...) =
    normalise!(deepcopy(t), args...; kwargs...)

"""
    remove_mean!(t::Trace) -> t

Remove the mean of trace `t` in place and return the trace.
"""
function remove_mean!(t::AbstractTrace)
    t.t .= t.t .- mean(t.t)
    t
end
remove_mean(t::AbstractTrace, args...; kwargs...) = remove_mean!(deepcopy(t), args...; kwargs...)

"""
    remove_trend!(t::Trace) -> t

Remove the trend from `t` in place and return the trace.
"""
function remove_trend!(t::AbstractTrace)
    time = times(t)
    x0, x1 = linear_regression(time, t.t)
    t.t .= t.t .- (x0 .+ x1.*time)
    t
end
remove_trend(t::AbstractTrace, args...; kwargs...) = remove_trend!(deepcopy(t), args...; kwargs...)
