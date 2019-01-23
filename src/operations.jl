# Trace operations

"""
    cut!(t, start, end; warn=true) -> t

Cut a `Trace` `t` in place between `start` and `end`.  An error is thrown if either
`start` or `end` are `missing`.

By default, a warning is shown if cut times lie outside the trace; set `warn` to `false`
to turn this off.
"""
function cut!(t::AbstractTrace, b, e; warn=true)
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
cut!(t::AbstractTrace, b::DateTime, e::DateTime; kwargs...) =
    cut(t, Dates.value(Dates.Nanosecond(b - t.evt.time))/1e9,
        Dates.value(Dates.Nanosecond(e - t.evt.time))/1e9; kwargs...)

"""
    cut!(t, pick1, offset1, pick2, offset; kwargs...) -> t
    cut!(t, pick, offset1, offset2; kwargs...) ->

Cut a trace `t` in place between `offset1` s after the first pick with name `pick1`
and `offset2` s after `pick2`.

In the second form, both offsets are relative to `pick`
"""
cut!(t::AbstractTrace, pick1, offset1, pick2, offset2; kwargs...) =
    cut!(t, first(picks(t, pick1)).time + offset1, first(picks(t, pick2)).time + offset2; kwargs...)
cut!(t::AbstractTrace, pick1, offset1, offset2; kwargs...) =
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
cut(t::AbstractTrace, args...; kwargs...) = cut!(deepcopy(t), args...; kwargs...)

# TODO: Enable antialiasing
"""
    decimate!(t, n; antialias=true) -> t
    decimate(t, n; antialias=true) -> t′

Decimate the trace `t` by removing all except every `n` points.  The sampling interval
is increased `n` times.  In the first form, update the trace in place
and return it.  In the second form, return an updated copy.

If `antialias` is `false`, then no antialiasing filtering is applied during decimation.
This means the decimated trace may contain spurious signals.
"""
function decimate!(t::AbstractTrace, n::Integer; antialias=true)
    antialias && error("Antialias filtering not implemented; use `antialias=false`")
    1 <= n || throw(ArgumentError("n must be grater than 0 (supplied $n)"))
    n == 1 && return t
    t.t = t.t[1:n:end]
    t.delta *= n
    t
end
decimate(t::AbstractTrace, n; antialias=true) = decimate!(deepcopy(t), n; antialias=antialias)
@doc (@doc decimate!) decimate

"""
    differentiate!(t::Trace; points=2) -> t
    differentiate(t::Trace; points=2) -> t′

Differentiate the trace `t` by performing `points`-point finite differencing.
In the first form, update the trace in place and return it.  In the second form,
return an updated copy.

### Available algorithms

- npoints == 2`: Two-point.  `dsdt.t[i] = (t.t[i+1] - t.t[i])/t.delta`.
  Non-central difference, so `t.b` is increased by half `t.delta`.
  The trace length is reduced by 1 samples.
- `points == 3`: Three-point. `dsdt.t[i] = (t.t[i+1] - t.t[i-1])/(2 * t.delta)`.
  Central difference.  `t.b` is increased by `t.delta`; the trace length is reduced
  by 2 samples.
- `points == 5`: Five-point. `dsdt.t[i] =
  (2/3)*(t.t[i+1] - t.t[i-1])/t.delta - (1/12)*(t.t[i+2] - t.t[i-2])/t.delta`.
  Central difference.  `t.b` is increased by `2t.delta`; `npts` reduced by 4.
"""
function differentiate!(t::AbstractTrace; points=2)
    points in (2, 3, 5) ||
        throw(ArgumentError("`points` must be one of (2, 3, 5)"))
    npts = nsamples(t)
    if points == 2
        @inbounds for i in 1:(npts-1)
            t.t[i] = (t.t[i+1] - t.t[i])/t.delta
        end
        pop!(t.t)
        t.b += t.delta/2
    elseif points == 3
        @inbounds for i in 2:(npts-1)
            t.t[i-1] = (t.t[i+1] - t.t[i-1])/(2*t.delta)
        end
        pop!(t.t); pop!(t.t)
        t.b += t.delta
    elseif points == 5
        t1 = (t.t[3] - t.t[1])/(2*t.delta)
        t2 = (t.t[end] - t.t[end-2])/(2*t.delta)
        d1 = 2/(3*t.delta)
        d2 = 1/(12*t.delta)
        t_minus_2 = t.t[1]
        t_minus_1 = t.t[2]
        tt = t.t[3]
        t_plus_1 = t.t[4]
        @inbounds for i in 2:(npts-3)
            t_plus_2 = t.t[i+3]
            t.t[i] = d1*(t_plus_1 - t_minus_1) - d2*(t_plus_2 - t_minus_2)
            t_minus_2 = t_minus_1
            t_minus_1 = tt
            tt = t_plus_1
            t_plus_1 = t_plus_2
        end
        t.t[1] = t1
        t.t[end-2] = t2
        pop!(t.t); pop!(t.t)
        t.b += t.delta
    end
    t
end
differentiate(t::AbstractTrace; kwargs...) = differentiate!(deepcopy(t); kwargs...)
@doc (@doc differentiate!) differentiate

"""
    envelope!(t::Trace) -> t
    envelope(t::Trace) -> t′

Replace the trace `t` with its envelope.
In the first form, update the trace in place and return the trace.
In the second form, return an updated copy.
"""
function envelope!(t::AbstractTrace)
    trace(t) .= abs.(DSP.hilbert(trace(t)))
    t
end
envelope(t::AbstractTrace) = envelope!(deepcopy(t))
@doc (@doc envelope!) envelope

"""
    flip!(t) -> t
    flip(t) -> t′

Reverse the direction of a trace so that it points the opposite way.

In the first form, update the trace in place and return the trace.
In the second form, return an updated copy.
"""
function flip!(t::AbstractTrace)
    any(ismissing, (t.sta.azi, t.sta.inc)) &&
        throw(ArgumentError("trace must have sta.azi and sta.inc defined"))
    t.sta.azi = mod(t.sta.azi+ 180, 360)
    t.sta.inc = 180 - t.sta.inc
    t.sta.cha = string(round(t.sta.azi, digits=2, base=10))
    trace(t)[:] .*= -1
    t
end
flip(t::AbstractTrace) = flip!(deepcopy(t))
@doc (@doc flip!) flip

"""
    integrate!(t::Trace, method=:trapezium) -> t
    integrate(t::Trace, method=:trapezium) -> t′

Replace `t` with its time-integral.  This is done by default using the trapezium rule.
Use `method=:rectangle` to use the rectangle rule.

In the first form, update the trace in place and return the trace.
In the second form, return an updated copy.

If `method==:trapezium` (the default), then the number of smples is reduced by one and
the begin time is increased by half the sampling interval.
"""
function integrate!(t::AbstractTrace, method::Symbol=:trapezium)
    npts = nsamples(t)
    if method == :trapezium
        total = zero(t.t[1])
        h = t.delta/2
        @inbounds for i in 1:(npts-1)
            total += h*(t.t[i] + t.t[i+1])
            t.t[i] = total
        end
        pop!(t.t)
        t.b += t.delta/2
    elseif method == :rectangle
        h = t.delta
        @inbounds for i in 2:npts
            t.t[i] = h*t.t[i] + t.t[i-1]
        end
    else
        throw(ArgumentError("`method` must by one of `:trapezium` or `:rectangle`"))
    end
    t
end
integrate(t::AbstractTrace, args...) = integrate!(deepcopy(t), args...)
@doc (@doc integrate!) integrate

"""
    normalise!(t::Trace, val=1) -> t
    normalise(t::Trace, val=1) -> t′

Normalise the trace `t` so that its maximum absolute amplitude is `val`.
In the first form, update the trace in place and return the trace.
In the second form, return an updated copy.

This function can also be spelled `normalize[!]`.
"""
function normalise!(t::AbstractTrace, val=1)
    maxval = maximum(abs, trace(t))
    t.t .*= val/maxval
    t
end
normalise(t::AbstractTrace, args...; kwargs...) =
    normalise!(deepcopy(t), args...; kwargs...)
@doc (@doc normalise!) normalise

LinearAlgebra.normalize!(t::AbstractTrace, args...; kwargs...) =
    normalise!(t, args..., kwargs...)
@doc (@doc normalise!) normalize!
LinearAlgebra.normalize(t::AbstractTrace, args...; kwargs...) =
    normalise(t, args...; kwargs...)
@doc (@doc normalise) = normalize

"""
    remove_mean!(t::Trace) -> t
    remove_mean(t::Trace) -> t′

Remove the mean of trace `t`.  In the first form, update the trace in place
and return it.  In the second form, return an updated copy.
"""
function remove_mean!(t::AbstractTrace)
    t.t .= t.t .- mean(t.t)
    t
end
remove_mean(t::AbstractTrace, args...; kwargs...) = remove_mean!(deepcopy(t), args...; kwargs...)
@doc (@doc remove_mean!) remove_mean

"""
    remove_trend!(t::Trace) -> t
    remove_mean(t::Trace) -> t′

Remove the trend from `t`.  In the first form, update the trace in place
and return it.  In the second form, return an updated copy.
"""
function remove_trend!(t::AbstractTrace)
    time = times(t)
    x0, x1 = linear_regression(time, t.t)
    t.t .= t.t .- (x0 .+ x1.*time)
    t
end
remove_trend(t::AbstractTrace, args...; kwargs...) = remove_trend!(deepcopy(t), args...; kwargs...)
@doc (@doc remove_trend!) remove_trend

"""
    taper!(t::AbstractTrace, width=0.05, form=:hanning) -> t
    taper(t::AbstractTrace, width=0.05, form=:hamming) -> t′

Apply a symmetric taper to each end of the data in trace `t`.
`form` may be one of `:hanning`, `:hamming` or `:cosine`.
`width` represents the fraction (at both ends) of the trace tapered, up to 0.5.

In the first form, update the trace in place and return it.  In the second form,
return an updated copy.
"""
function taper!(t::AbstractTrace, width=0.05; form::Symbol=:hanning)
    form in (:hamming, :hanning, :cosine) ||
        throw(ArgumentError("`form` must be one of `:hamming`, `:hanning` or `:cosine`"))
    0 < width <= 0.5 || throw(ArgumentError("SAC.taper!: width must be between 0 and 0.5"))
    n = max(2, floor(Int, (nsamples(t) + 1)*width))

    T = eltype(trace(t))
    npts = nsamples(t)

    if form in (:hamming, :hanning)
        omega = T(π/n)
        if form == :hanning
            f0 = f1 = T(0.50)
        elseif form == :hamming
            f0 = T(0.54)
            f1 = T(0.46)
        end

        @inbounds for i in 0:n-1
            amp = f0 - f1*cos(omega*T(i))
            j = npts - i
            t.t[i+1] *= amp
            t.t[j] *= amp
        end
    end

    if form == :cosine
        omega = T(π/2n)
        @inbounds for i in 0:n-1
            amp = sin(omega*i)
            j = npts - i
            t.t[i+1] *= amp
            t.t[j] *= amp
        end
    end

    t
end
taper(t::AbstractTrace, args...; kwargs...) = taper!(deepcopy(t), args...; kwargs...)
@doc (@doc taper!) taper
