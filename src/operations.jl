# Trace operations

"""
    cut!(t, start, end; allowempty=false, warn=true) -> t

Cut a `Trace` `t` in place between `start` and `end`.  An error is thrown if either
`start` or `end` are `missing`.

An error is thrown if the trace would be empty because either the end cut time is
before the start of the trace, or the start cut is after the end, unless `allowempty`
is `true`.

By default, a warning is shown if cut times lie outside the trace; set `warn` to `false`
to turn this off.

# Example
```
julia> t = Trace(0, 1, [0, 1, 2, 3, 4, 5]);

julia> trace(cut!(t, 2, 4))
3-element Array{Float64,1}:
 2.0
 3.0
 4.0
```
---
    cut!(t, start_date, end_date; kwargs...) -> t

Cut a `Trace` `t` in place between dates `start_date` and `end_date`.

# Example
Create a trace starting at midnight on 1 January 3000, and cut to
between 10 s and 1 minute after midnight:
```
julia> using Dates: DateTime

julia> t = Trace(0, 1, 100); # 1 Hz sampling for 100 s;

julia> t.evt.time = DateTime(3000, 1, 1); # The year 3000;

julia> cut!(t, DateTime(3000, 1, 1, 0, 0, 10), DateTime(3000, 1, 1, 0, 1, 0))
Seis.Trace{Float64,Vector{Float64},Seis.Geographic{Float64}}:
            b: 10.0
        delta: 1.0
 GeogStation{Float64}:
     sta.meta: Seis.SeisDict{Symbol, Any}()
 GeogEvent{Float64}:
     evt.time: 3000-01-01T00:00:00
     evt.meta: Seis.SeisDict{Symbol, Any}()
 Trace:
        picks: 0
         meta: 

julia> startdate(t), enddate(t)
(DateTime("3000-01-01T00:00:10"), DateTime("3000-01-01T00:01:00"))
```
"""
function cut!(t::AbstractTrace, b, e; allowempty=false, warn=true)
    ib, ie, b′, empty = _cut_time_indices(t, b, e; allowempty=allowempty, warn=warn)
    if empty
        empty!(trace(t))
    else
        t.t = t.t[ib:ie]
    end
    t.b = b′
    t
end

"""
    cut!(t, pick1, offset1, pick2, offset; kwargs...) -> t
    cut!(t, pick, offset1, offset2; kwargs...) ->

Cut a trace `t` in place between `offset1` s after the first pick `pick1`
and `offset2` s after `pick2`.

In the second form, both offsets are relative to `pick`.

The values of `pick1`, `pick2` and `pick` are passed to [`picks`](@ref)
and so may be a `Symbol` (giving the key of the pick), a `String` (giving
the pick name) or a `Regex` (which matches the pick name).

# Example
```
julia> t = sample_data();

julia> starttime(t), endtime(t)
(52.66f0, 62.65f0)

julia> cut!(t, :A, 0, :F, 1);

julia> starttime(t), endtime(t)
(53.67f0, 61.979996f0)
```
"""
function cut!(t::AbstractTrace, pick1, offset1, pick2, offset2; kwargs...)
    picks1 = picks(t, pick1, sort=:time)
    picks2 = picks(t, pick2, sort=:time)
    isempty(picks1) && throw(ArgumentError("trace does not contain pick '$pick1'"))
    isempty(picks2) && throw(ArgumentError("trace does not contain pick '$pick2'"))
    cut!(t, first(picks1).time + offset1,
            first(picks2).time + offset2; kwargs...)
end

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

See also: [`cut!`](@ref), [`picks`](@ref).
"""
cut(t::AbstractTrace, args...; kwargs...) = cut!(deepcopy(t), args...; kwargs...)

"""
    _cut_time_indices(t::AbstractTraceArray, b, e; warn, allowempty) -> ib, ie, b′, empty

Compute the start index `ib` and end index `ie` which cuts the trace array
`t` in time between `b` and `e`.  If the trace is empty, `empty` is `false`.
It is otherwise `true`.  Also return the adjusted trace start time, b′.

`b` and `e` may be both either time in s relative to the trace array
event time (if any), or absolute `DateTime`s.  In the latter case,
`t.evt.time` must be set.

!!! note
    Users implementing methods of [`cut`](@ref) and [`cut!`](@ref)
    can make use of this function, but note that it is subject to
    change.  [`nearest_sample`](@ref) can be used to implement the
    logic if absolute stability is required.
"""
function _cut_time_indices(t::AbstractTrace, b, e; warn=true, allowempty=false)
    (b === missing || e === missing) && throw(ArgumentError("Start or end cut time is `missing`"))
    e < b && throw(ArgumentError("End cut time ($e) is before start cut ($b)"))

    # Start and end times of trace before cutting
    old_start = starttime(t)
    old_end = endtime(t)

    # First and last samples of data
    istart = first(first(axes(trace(t))))
    iend = last(first(axes(trace(t))))

    # Samples at which to cut
    ib = round(Int, (b - old_start)/t.delta) + 1
    ib = max(ib, istart)
    ie = iend - round(Int, (old_end - e)/t.delta)
    ie = min(ie, iend)

    new_start = old_start + (ib - istart)*t.delta

    if b > endtime(t) || e < starttime(t)
        empty = true
        if !allowempty
            b > endtime(t) &&
                throw(ArgumentError("Beginning cut time $b is later than end of trace ($old_end)."))
            e < starttime(t) &&
                throw(ArgumentError("End cut time $e is earlier than start of trace ($old_start)."))
        end
        new_start = b
    else
        empty = false
        if b < starttime(t)
            warn && @warn("Beginning cut time $b is before start of trace.  Setting to $old_start.")
            new_start = old_start
        end
        if e > endtime(t)
            warn && @warn("End cut time $e is after end of trace.  Setting to $old_end.")
        end
    end

    ib, ie, new_start, empty
end

function _cut_time_indices(t::AbstractTrace, b::Dates.AbstractDateTime, e::Dates.AbstractDateTime; warn=true, allowempty=true)
    t.evt.time === missing && throw(ArgumentError("no event time defined for trace array"))
    b_ns::Dates.Nanosecond = b - origin_time(t)
    e_ns::Dates.Nanosecond = e - origin_time(t)
    b = Dates.value(b_ns)/1_000_000_000
    e = Dates.value(e_ns)/1_000_000_000
    _cut_time_indices(t, b, e; warn, allowempty)
end

"""
    decimate!(t, n; antialias=true) -> t
    decimate(t, n; antialias=true) -> t′

Decimate the trace `t` by removing all except every `n` points.  The sampling interval
is increased `n` times.  In the first form, update the trace in place
and return it.  In the second form, return an updated copy.

By default, an antialiasing and decimation FIR filter is applied.  This may cause
artifacts in the signal at the extremes of the trace.

If `antialias` is `false`, then no antialiasing filtering is applied during decimation.
This means the decimated trace may contain spurious signals.
"""
function decimate!(t::AbstractTrace, n::Integer; antialias=true)
    1 <= n || throw(ArgumentError("n must be greater than 0 (supplied $n)"))
    n == 1 && return t
    if antialias
        t.t = DSP.resample(t.t, 1//n)
    else
        t.t = t.t[1:n:end]
    end
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

- `points == 2`: Two-point.  `dsdt.t[i] = (t.t[i+1] - t.t[i])/t.delta`.
  Non-central difference, so `t.b` is increased by half `t.delta`.
  The trace length is reduced by 1 samples.
- `points == 3`: Three-point. `dsdt.t[i] = (t.t[i+1] - t.t[i-1])/(2 * t.delta)`.
  Central difference.  `t.b` is increased by `t.delta`; the trace length is reduced
  by 2 samples.
- `points == 5`: Five-point. `dsdt.t[i] =
  (2/3)*(t.t[i+1] - t.t[i-1])/t.delta - (1/12)*(t.t[i+2] - t.t[i-2])/t.delta`,
  except for the first and last points, which use a three-point central difference
  meaning only two points fewer are retained as for `points == 3`.
  Central difference.  `t.b` is increased by `t.delta`; `npts` reduced by 2.

# Example
```
julia> t = Trace(0, 1, [0, 1, -1, 0]);

julia> d = differentiate(t); trace(d)
3-element Array{Float64,1}:
  1.0
 -2.0
  1.0

julia> starttime(d)
0.5
```
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

# Example
```
julia> t = Trace(0, 1, [0, 0, 0, 1, -1, 0, 0, 0]);

julia> trace(envelope(t))
8-element Array{Float64,1}:
 0.10355339059327379
 0.10355339059327379
 0.6035533905932737
 1.1680225577002512
 1.1680225577002512
 0.6035533905932737
 0.10355339059327373
 0.10355339059327379
```
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
This preserves the sense of the data; for example, a positive signal on
an eastward-pointing channel becomes a negative signal on the flipped
westward pointing channel.  Both before and after, the signal is
positive eastwards.

The `t.sta` must contain both azimuth and inclination information.

In the first form, update the trace in place and return the trace.
In the second form, return an updated copy.

# Example
```
julia> t = Trace(0, 1, [0, 1, 0]); # Positive arrival at 1 s

julia> t.sta.azi, t.sta.inc = 0, 90 # North horizontal component
(0, 90)

julia> flip!(t)
Seis.Trace{Float64,Array{Float64,1},Seis.Geographic{Float64}}:
            b: 0.0
        delta: 1.0
 Station{Float64,Seis.Geographic{Float64}}:
      sta.cha: 180.0
      sta.azi: 180.0
      sta.inc: 90.0
     sta.meta: Seis.SeisDict{Symbol,Any}()
 Event{Float64,Seis.Geographic{Float64}}:
     evt.meta: Seis.SeisDict{Symbol,Any}()
 Trace:
        picks: 0
         meta: 

julia> trace(t)
3-element Array{Float64,1}:
 -0.0
 -1.0
 -0.0
```
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

If `method==:trapezium` (the default), then the number of samples is reduced by one and
the begin time is increased by half the sampling interval.

# Example
```
julia> t = Trace(0, 0.1, [0, 1, 1, 0]);

julia> trace(integrate(t))
3-element Array{Float64,1}:
 0.05
 0.15000000000000002
 0.2

julia> trace(integrate(t, :rectangle))
4-element Array{Float64,1}:
 0.0
 0.1
 0.2
 0.2
```
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
    merge!(t1::AbstractTrace, ts::AbstractArray{<:AbstractTrace}; gaps=:zero, overlaps=:mean, sample_tol=0.1, check=true) -> t1
    merge!(t1, ts...; kwargs...) -> t1
    merge!([t1, ts...]; kwargs...) -> t1

    merge(t1::AbstractTrace, ts::AbstractArray{<:AbstractTrace}; gaps=:zero, overlaps=:mean, sample_tol=0.1, check=true) -> t1
    merge(t1, ts...; kwargs...) -> t1
    merge([t1, ts...]; kwargs...) -> t1

Merge two or more traces together into the first trace, retaining only
station and event information from the first trace.

In the first form, update the trace in place and return the trace.
In the second form, return an updated copy.

All traces must have the same sampling interval.  They must also have the
same channel code (see [`channel_code`](@ref)), unless `check` is `false`.
The traces do not need to have the same element type or geometry.

Where gaps between traces occur, these may be filled in a number of
different ways, or an error may be thrown.  Tapering may be applied to
the data either side of gaps.  See 'Keyword arguments' below for details.

Likewise, where overlaps occur, either data from the first or last trace
in each overlap may be used, or the mean of the traces used.  Other
options are also possible.

If the start times of traces (and hence each sample itself) are not
quantised in time in the same way, then an error is thrown, unless
the difference in quantisation is less than a fraction of `sample_tol`
of the sampling interval.  Hence `sample_tol` should be between 0 and
0.5.

If all traces have an event time set, then merging is done in absolute
time.  If not all have an event time set, then all times are assumed
to be relative to the same origin and only the traces' `b` fields are
used.  To perform merging in relative time regardless of whether the
origin time is set, pass `relative = true`.

Empty traces are checked for matching channel codes, but are not otherwise
merged into the final trace.  If all traces are empty, the first is returned
unaltered.

# Keyword arguments

## `gaps`
`gaps` controls how gaps between continuous segments of data are
treated, and may take one of the following values:
- `:zero` (default): Fill any gaps with zero
- `:error`: Throw an error if any gaps are present
- `:linear`: Linearly interpolate between the last sample before the gap
  and the first sample after the gap.  Cannot be used with tapering.
- `value`: Fill with `value`, which must be convertable to the element
  type of `t1`.

## `overlaps`
`overlaps` controls how the new merged trace uses traces which overlap
in time, and may take one of the following values:
- `:mean` (default): Take the mean of the values at each overlapping
  sample
- `:first`: Use the data from the first (in time) trace in the overlap
- `:last`: Use the data from the last (in time) trace
- `:zero`: Zero out any overlapping periods
- `:error`: Throw an error if any overlaps are present and the data
  are not identical.  (Note that the error type is not at present defined
  but may be in a future version.)
- `value`: Fill with `value`, which must be convertable to the element
  type of `t1`.

## `sample_tol`
Any trace which does not have its samples quantised to the same
as `t1` to within a fraction of a sampling interval `sample_tol`
will cause an error to be thrown.  Set this to a value between 0 and 0.5
to control this.  A value of 0.5 will mean any quantisation differences
are ignored.

## `taper`
If `taper` is set to a fractional width, a taper is performed
on the ends of traces adjoining gaps, with `taper` defining the
proportional length of the taper relative to the gap length.  Values are
tapered to 0.

Note that if `taper` is used with a number for `gaps`, then sharp
jumps will occur at the first and last sample of each gap to whatever
`value` is supplied to the `gaps` argument

## `taper_form`
Determines the type of taper applied around gaps if `taper` is set.
This can be one of `:hanning` (the default), `:hamming` or `:cosine`.
See [`taper!`](@ref) for details.

## `relative`
If `relative` is `true`, then ignore any origin times in traces
and merge them all only with regard to their [`starttime`](@ref).
By default merging is done in absolute time if all traces have `.evt.time`
set.

## `metadata`
If `true` (the default), then entries in the `.meta` field of each trace
are merged into the `.meta` field of the first trace.  Entries of the first
trace are preserved, while entries not in the first trace are added in turn
from the last trace to the first.  This means duplicate entries are not
overwritten in the first trace, and entries present in more than one of
the other traces are taken from the second, third, etc., trace in preference
to any later traces.
"""
function Base.merge!(
    t1::AbstractTrace,
    ts::Union{AbstractArray{<:AbstractTrace}, NTuple{N,AbstractTrace} where N};
    gaps=:zero,
    overlaps=:mean,
    taper=0,
    taper_form=:hanning,
    relative=false, check=true,
    metadata=true,
    sample_tol=0.1
)
#=
There are several possibilities in which two traces may line up:

Perfect overlap:
    t1: |--------------------|
    t2: |--------------------|
                                      --------> time
Complete gap
    t1: |---------------|
    t2:                      |--------------------|

    t1:                      |--------------------|
    t2: |---------------|

Partial overlap
    t1:           |--------------------|
    t2:                   |--------------------|

    t1:                   |--------------------|
    t2:           |--------------------|
 
Containment
    t1:           |--------------------|
    t2:                   |------|

    t1:                   |------|
    t2:           |--------------------|

Our method simply takes the earliest and latest times, sorts the traces by
start time, creates a new trace covering the whole time span, and fills in
the new trace, looking for places where the above cases need to be handled:

    t1: |--------------------|
    t2:                          |----|
    t3:                             |-----|

        |<- start time         end time ->|
                     samples
        |||||||||||||||||||||||||||||||||||

                   merged trace
        1111111111111111111111___222***3333
               (_: gap; *: overlap)
=#
    0 <= sample_tol <= 0.5 ||
        throw(ArgumentError("sample_tol should be between 0 and 0.5"))

    all(t -> t.delta == t1.delta, ts) ||
        throw(ArgumentError("all traces must have the same sampling interval"))
    delta = t1.delta

    if check
        _channel_codes_are_equal(t1, ts) ||
            throw(ArgumentError("all traces must have the same channel code"))
    end

    ts_nonempty = filter(!iszero∘nsamples, ts)

    # Cannot both linearly interpolate and taper
    if gaps === :linear && !iszero(taper)
        throw(ArgumentError("cannot combine `gaps=:linear` and `taper>0`"))
    end

    0 <= taper ||
        throw(ArgumentError("taper must be positive"))
    taper_form in (:hanning, :hamming, :cosine) ||
        throw(ArgumentError("taper_form must be one of :hanning, :hamming or :cosine"))

    # Do it by date unless any of the traces doesn't have absolute time,
    # in which case do it all by relative time, or if we have asked for
    # relative time
    time_base = (relative ||
            (t1.evt.time === missing || any(t -> t.evt.time === missing, ts_nonempty))) ?
        :relative : :absolute

    # These update and return `t1`
    if time_base == :absolute
        _merge_absolute!(t1, ts_nonempty, gaps, overlaps, taper, sample_tol, taper_form)
    elseif time_base == :relative
        _merge_relative!(t1, ts_nonempty, gaps, overlaps, taper, sample_tol, taper_form)
    end

    # Merge the metadata
    if metadata
        for t in ts
            t1.meta = merge(t.meta, t1.meta)
        end
    end

    t1
end

function Base.merge!(t1::AbstractTrace, t2::AbstractTrace, ts::Vararg{AbstractTrace}; kwargs...)
    merge!(t1, [t2, ts...]; kwargs...)
end
function Base.merge!(ts::AbstractArray{<:AbstractTrace}; kwargs...)
    isempty(ts) && throw(ArgumentError("cannot merge an empty array of traces"))
    length(ts) == 1 && return only(ts)
    merge!(ts[begin], @view(ts[begin+1:end]); kwargs...)
end

"""
    merge(t1::AbstractTrace, t2::AbstractTrace...) -> t_merged
    merge(t1::AbstractTrace, ::Vararg{AbstractTrace}) -> t_merged
    merge(ts::AbstractArray{<:AbstractTrace}) -> t_merged

For details of out-of-place trace merging, see [`merge!`](@ref merge!(::AbstractTrace, ::AbstractArray{<:AbstractTrace})).
"""
Base.merge(t1::AbstractTrace, ts::AbstractArray{<:AbstractTrace}; kwargs...) =
    merge!(deepcopy(t1), ts; kwargs...)
Base.merge(t1::AbstractTrace, t2::AbstractTrace, ts::Vararg{AbstractTrace}; kwargs...) =
    merge!(deepcopy(t1), [t2, ts...]; kwargs...)
Base.merge(t1::AbstractTrace, ts::NTuple{N,AbstractTrace} where N; kwargs...) =
    merge!(deepcopy(t1), ts...; kwargs...)
function Base.merge(ts::AbstractArray{<:AbstractTrace}; kwargs...)
    isempty(ts) && throw(ArgumentError("cannot merge an empty array of traces"))
    length(ts) == 1 && return deepcopy(only(ts))
    merge!(deepcopy(ts[begin]), @view(ts[begin+1:end]); kwargs...)
end

"""
    _merge_relative!(t1, ts, gaps, overlaps, taper, sample_tol, taper_form) -> t1

Merge all the data in the set of traces `ts` into the single `Trace` `t1`.

See `Base.merge!(::AbstractTrace, ::AbstractArray{<:AbstractTrace})` for
details of the arguments `gaps`, `overlaps`, `taper`, `sample_tol` and `taper_form`.
"""
function _merge_relative!(t1, ts, gaps, overlaps, taper, sample_tol, taper_form)
    # Just put points in t1 if there is only one non-empty trace at `t1` is empty
    if nsamples(t1) == 0 && length(ts) == 1
        t_nonempty = only(ts)
        append!(trace(t1), trace(t_nonempty))
        t1.b = t_nonempty.b
        return t1
    elseif nsamples(t1) == 0 && isempty(ts)
        return t1
    end

    # Already checked they all have same delta
    delta = t1.delta

    # Check for quantisation of samples
    t1_offset = _round_offset(starttime(t1), delta)

    if any(t -> abs(_round_offset(starttime(t), delta)) > sample_tol*delta, ts)
        throw(ArgumentError(
            "trace sample times are not aligned to within $(sample_tol) of a sample ($(sample_tol*delta) s)"))
    end

    # Remove the first trace if it's empty
    all_ts = nsamples(t1) == 0 ? ts : [t1; ts]

    # Sort data
    sorted_inds = sortperm(all_ts, by=starttime)

    # New start time for all traces based on t1, but quantised the same as t1
    first_sample = _quantise(minimum(starttime, all_ts), delta, t1_offset)
    # Final sample of new data, quantised same as t1
    last_sample = _quantise(maximum(endtime, all_ts), delta, t1_offset)

    # Simple case of all traces ordered without gaps or overlaps
    if _traces_are_sequential(all_ts, sorted_inds, sample_tol)
        t1.t = reduce(vcat, trace(all_ts[i]) for i in sorted_inds)
        t1.b = first_sample
        return t1
    end

    # New raw data for merged `t1`
    new_nsamples = round(Int, (last_sample - first_sample)/delta) + 1
    new_axes = (new_nsamples, Base.tail(axes(trace(t1)))...)
    new_trace = similar(trace(t1), new_axes)
    new_trace .= 0

    # Trace keeping track of number of original traces at each sample
    count_trace = zeros(Int, new_nsamples)

    # Fill data into new trace
    for i in sorted_inds
        t = all_ts[i]
        isempty(trace(t)) && continue

        b = _quantise(starttime(t), delta, t1_offset)
        i1 = round(Int, (b - first_sample)/delta) + 1
        i2 = i1 + nsamples(t) - 1

        if !(i1 in eachindex(new_trace) || i2 in eachindex(new_trace))
            error("error in indexing logic; please report bug")
        end

        if overlaps === :mean
            @views new_trace[i1:i2] .+= trace(t)
        elseif overlaps === :first
            @views new_trace[i1:i2] .= ifelse.(count_trace[i1:i2] .> 0, new_trace[i1:i2], trace(t))
            # for (i, new_i) in enumerate(i1:i2)
            #     new_trace[new_i] = count_trace[new_i] > 0 ? new_trace[new_i] : trace(t)[begin+i-1]
            # end
        elseif overlaps === :last
            new_trace[i1:i2] .= trace(t)
        elseif overlaps === :error
            any(i -> count_trace[i] > 1, i1:i2) &&
                error("traces overlap")
        elseif !isa(overlaps, Symbol) # For zero and user-supplied value cases, sort it out later
            new_trace[i1:i2] .= trace(t)
        else
            overlaps in (:zero,) ||
                throw(ArgumentError("unrecognised overlaps option '$overlaps'"))
        end

        @views count_trace[i1:i2] .+= 1
    end

    has_overlaps = any(>(1), count_trace)

    # Sort out overlaps
    if has_overlaps
        for i in eachindex(new_trace, count_trace)
            if count_trace[i] > 1
                if overlaps === :mean
                    new_trace[i] /= max(count_trace[i], 1)
                elseif overlaps === :zero
                    new_trace[i] = count_trace[i] > 1 ? zero(eltype(new_trace)) : new_trace[i]
                elseif !isa(overlaps, Symbol) # User-supplied value
                    new_trace[i] = count_trace[i] > 1 ? overlaps : new_trace[i]
                end
            end
        end
    end

    # Sort out gaps and apply tapers (`:zero` case already handled)
    gap_sections = _find_gap_sections(count_trace)
    if !isempty(gap_sections)
        if gaps === :error
            error("gaps present in merged trace")
        elseif gaps === :zero
            # Case already handled, but keep this here to catch
            # unsupported options at the end of this `if` block
        elseif gaps === :linear
            for (i1, i2) in gap_sections
                if i1 <= firstindex(new_trace) || i2 >= lastindex(new_trace)
                    error("logic error in gap filling code; please report bug")
                end
                # Single sample gap
                if i1 == i2
                    new_trace[i1] = (new_trace[i1-1] + new_trace[i2+1])/2
                else
                    v1 = new_trace[i1-1]
                    v2 = new_trace[i2+1]
                    nsamples_in_gap = i2 - i1 + 1
                    new_trace[i1-1:i2+1] .= range(v1, v2, length=nsamples_in_gap+2)
                end
            end
        elseif !isa(gaps, Symbol)
            for (i1, i2) in gap_sections
                new_trace[i1:i2] .= gaps
            end
        else
            throw(ArgumentError("unrecognised gaps option '$gaps'"))
        end

        if !iszero(taper)
            for (igap, (i1, i2)) in enumerate(gap_sections)
                gap_nsamples = i2 - i1 + 1
                taper_nsamples = round(Int, taper*gap_nsamples)

                if igap == 1
                    samples_before = firstindex(new_trace):(i1 - 1)
                    n_before = min(length(samples_before), taper_nsamples)
                else
                    samples_before = (gap_sections[igap-1][2] + 1):(i1 - 1)
                    n_before = min(length(samples_before)÷2, taper_nsamples)
                end

                _taper_core!(view(new_trace, samples_before), n_before, false, true, taper_form)

                if igap == lastindex(gap_sections)
                    samples_after = (i2 + 1):lastindex(new_trace)
                    n_after = min(length(samples_after), taper_nsamples)
                else
                    samples_after = (i2 + 1):(gap_sections[igap+1][1] - 1)
                    n_after = min(length(samples_after)÷2, taper_nsamples)
                end

                _taper_core!(view(new_trace, samples_after), n_after, true, false, taper_form)
            end
        end
    end

    t1.b = first_sample
    t1.t = new_trace

    t1
end

function _merge_absolute!(t1, ts, gaps, overlaps, taper, sample_tol, taper_form)
    # FIXME: Make this more efficient; maybe by passing the relative start times
    shifted_ts = origin_time.(ts, t1.evt.time)
    _merge_relative!(t1, shifted_ts, gaps, overlaps, taper, sample_tol, taper_form)
end

"""
Return `true` if the channel codes are equal for both traces.
"""
function _channel_codes_are_equal(t1::AbstractTrace, t2::AbstractTrace)
    isequal(t1.sta.net, t2.sta.net) &&
        isequal(t1.sta.sta, t2.sta.sta) &&
        isequal(t1.sta.loc, t2.sta.loc) &&
        isequal(t1.sta.cha, t2.sta.cha)
end
_channel_codes_are_equal(t1::AbstractTrace, ts) =
    all(_channel_codes_are_equal(t1, t) for t in ts)

"""
Find an offset of a time relative to a timing base specified by samples
at `delta` s increments, starts at 0 s.
"""
function _round_offset(start_time, delta)
    offset = start_time%delta
    if offset < -delta/2
        offset + delta
    elseif offset > delta/2
        offset - delta
    else
        offset
    end
end

"""
For the set of traces `ts`, whose ordered starttimes are given by `indices`,
return `true` is each trace succeeds the previous one and the first sample
of each trace is one sample later than the last sample of the previous;
otherwise return `false`.
"""
function _traces_are_sequential(ts, indices, sample_tol)
    ntraces = length(ts)
    length(indices) == ntraces || error("`ts` and `indices` must be the same length")
    ntraces > 1 ||
        throw(ArgumentError("_traces_are_sequential should not get called with less than two traces"))

    delta = first(ts).delta
    delta_min, delta_max = (1 - sample_tol)*delta, (1 + sample_tol)*delta

    for index in eachindex(indices)
        index == lastindex(indices) && return true
        if delta_min <= (starttime(ts[indices[index + 1]])
                - endtime(ts[indices[index]])) <= delta_max
            continue
        else
            return false
        end
    end
end

"""
Given `count_trace`, a vector containing the number of traces which
have contributed to each sample in a final merged trace, return
a `Vector{Tuple{Int,Int}}`, where each element of the vector contains
the start and end index of `count_trace` where a run of `0`s occur.
In other words, each pair of indices gives the start and end of a gap
in the merged trace.
"""
function _find_gap_sections(count_trace)
    gap_start = typemin(Int) + 1
    gap_end = typemin(Int)
    gap_sections = Tuple{Int,Int}[]

    for (i, n) in pairs(count_trace)
        if n == 0
            if i == gap_end + 1
                # Gap carries on
                gap_end = i
            else
                # New gap
                gap_start = i
                gap_end = i
            end
        end

        if n > 0 || i == lastindex(count_trace)
            # Not a gap or the last sample
            if gap_end >= gap_start
                push!(gap_sections, (gap_start, gap_end))
                gap_start = typemin(Int) + 1
                gap_end = typemin(Int)
            end
        end
    end

    gap_sections
end

"""
    normalise!(t::Trace, val=1) -> t
    normalise(t::Trace, val=1) -> t′

Normalise the trace `t` so that its maximum absolute amplitude is `val`.
In the first form, update the trace in place and return the trace.
In the second form, return an updated copy.

This function can also be spelled `normalize[!]`.

# Example
```
julia> t = Trace(0, 0.1, [0, -1, 2]);

julia> trace(normalise(t))
3-element Array{Float64,1}:
  0.0
 -0.5
  1.0

julia> trace(normalise(t, 2))
3-element Array{Float64,1}:
  0.0
 -1.0
  2.0
```
"""
function normalise!(t::AbstractTrace, val=1)
    maxval = maximum(abs, trace(t))
    iszero(maxval) && return t
    t.t .*= val/maxval
    t
end
normalise(t::AbstractTrace, args...; kwargs...) =
    normalise!(deepcopy(t), args...; kwargs...)
@doc (@doc normalise!) normalise

# Avoid method ambiguity with LinearAlgebra.normalize(x, p::Real)
LinearAlgebra.normalize!(t::AbstractTrace, val::Real=1) =
    normalise!(t, val)
@doc (@doc normalise!) normalize!
LinearAlgebra.normalize(t::AbstractTrace, val::Real=1) =
    normalise(t, val)
@doc (@doc normalise) normalize

"""
    remove_mean!(t::Trace) -> t
    remove_mean(t::Trace) -> t′

Remove the mean of trace `t`.  In the first form, update the trace in place
and return it.  In the second form, return an updated copy.

# Example
```
julia> t = Trace(0, 0.01, [1, 1, 3, -1]);

julia> trace(remove_mean(t))
4-element Array{Float64,1}:
  0.0
  0.0
  2.0
 -2.0
```
"""
function remove_mean!(t::AbstractTrace)
    t.t .= t.t .- mean(t.t)
    t
end
remove_mean(t::AbstractTrace, args...; kwargs...) = remove_mean!(deepcopy(t), args...; kwargs...)
@doc (@doc remove_mean!) remove_mean

"""
    remove_trend!(t::Trace) -> t
    remove_trend(t::Trace) -> t′

Remove the trend from `t`.  In the first form, update the trace in place
and return it.  In the second form, return an updated copy.

# Example
```
julia> t = Trace(0, 0.2, [1, 2, 3, 4]);

julia> trace(remove_trend(t))
4-element Array{Float64,1}:
 -2.220446049250313e-16
  0.0
  0.0
  4.440892098500626e-16
```
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
    resample!(t::AbstractTrace; delta, n) -> t
    resample(t::AbstractTrace; delta, n) -> t′

Resample the trace `t` so that either the sampling interval becomes `delta` s,
or its sampling rate is increased `n` times.  One of `delta` or `n` must
be given.

In the first form, update the trace in place and return it.  In the second form,
return an updated copy (`t′`).

The functions uses [`DSP.resample`](@ref) to perform the operation, which applies
an antialias filter and 'additional operations' to prevent aliasing and minimise
other artifacts.

To perform decimation without antialiasing, use [`decimate`](@ref) or [`decimate!`](@ref)
with `antialias=false`.

# Example
```
julia> t = Trace(0, 0.5, 1:4);

julia> trace(resample(t, 3))
```
"""
function resample!(t::AbstractTrace; delta=nothing, n=nothing)
    if (delta === nothing && n === nothing) || (delta !== nothing && n !== nothing)
        throw(ArgumentError("one and only of `delta` and `n` must be given"))
    end
    rate = if delta !== nothing
        delta == t.delta && return t # Nothing to do
        t.delta/delta
    else # n !== nothing
        n == 1 && return t # Nothing to do
        n
    end
    newtrace = DSP.resample(trace(t), rate)
    t.t = newtrace
    t.delta = t.delta/rate
    t
end
DSP.resample(t::AbstractTrace; kwargs...) = resample!(deepcopy(t); kwargs...)
@doc (@doc resample!) resample


"""
    taper!(t::AbstractTrace, width=0.05; left=true, right=true, form=:hanning, time=nothing) -> t
    taper(t::AbstractTrace, width=0.05; left=true, right=true, form=:hamming, time=nothing) -> t′

Apply a taper to the ends of the data in trace `t`.
`form` may be one of `:hanning`, `:hamming` or `:cosine`.
`width` represents the fraction (at both ends) of the trace tapered, up to 0.5.

Optionally, specify `time` as an absolute length in time for the tapering
period (at both ends), in which case `width` is ignored.

By default, tapering is applied to both ends.  If `left` is `false`, then only
the 'right' end (later in time part) of the trace is tapered; and vice versa.

In the first form, update the trace in place and return it.  In the second form,
return an updated copy.

# Example
```
julia> t = Trace(0, 1, [-1, 1, -1, 1, -1, 1]);

julia> trace(taper(t))
6-element Array{Float64,1}:
 -0.0
  0.49999999999999994
 -1.0
  1.0
 -0.49999999999999994
  0.0

julia> trace(taper(t; time=0.25))
6-element Vector{Float64}:
 -0.0
  0.49999999999999994
 -1.0
  1.0
 -0.49999999999999994
  0.0

julia> trace(taper(t; right=false))
6-element Vector{Float64}:
 -0.0
  0.49999999999999994
 -1.0
  1.0
 -1.0
  1.0
```
"""
function taper!(t::AbstractTrace, width=0.05;
        left::Bool=true, right::Bool=true, form::Symbol=:hanning, time=nothing)
    form in (:hamming, :hanning, :cosine) ||
        throw(ArgumentError("`form` must be one of `:hamming`, `:hanning` or `:cosine`"))

    left == right == false && return t

    if time !== nothing
        iszero(time) && return t
        times = Seis.times(t)
        trace_length = last(times) - first(times)
        width = time/trace_length
        if left && right
            0 < width <= 0.5 ||
                throw(ArgumentError(
                    "time must be between 0 s and half the trace length for a two-sided taper"))
        else
            0 < width <= 1 ||
                throw(ArgumentError(
                    "time must be between 0 s and the trace length for one-sided taper"))
        end
    else
        iszero(width) && return t
        if left && right
            0 < width <= 0.5 ||
                throw(ArgumentError("width must be between 0 and 0.5 for a two-sided taper"))
        else
            0 < width <= 1 ||
                throw(ArgumentError("time must be between 0 s and the trace length for a one-sided taper"))
        end
    end

    n = max(
        2,
        min(
            floor(Int, (nsamples(t) + 1)*width),
            nsamples(t)
        )
    )

    _taper_core!(trace(t), n, left, right, form)

    t
end
taper(t::AbstractTrace, args...; kwargs...) = taper!(deepcopy(t), args...; kwargs...)
@doc (@doc taper!) taper

"""
    _taper_core!(data, n, left, right, form)

Perform tapering using a window of kind `form`, on `n` samples of `data`.  The
`left` and `right` are `Bool`s saying whether or not respectively to taper the
first and last part of `data`.

Note that this works when `data` is a `AbstractVector`, or `AbstractMatrix` where
the individual traces are arranged down columns.
"""
function _taper_core!(data, n, left, right, form)
    T = eltype(data)

    n <= size(data, 1) || throw(BoundsError(data, n))

    if form in (:hamming, :hanning)
        # 1/n, without the π factor (since we call cospi)
        omega_π = 1/T(n)
        if form == :hanning
            f0 = f1 = T(0.50)
        elseif form == :hamming
            f0 = T(0.54)
            f1 = T(0.46)
        end

        for j in axes(data, 2)
            @inbounds for i in 0:(n - 1)
                amp = f0 - f1*cospi(omega_π*T(i))
                data[begin+i,j] *= left ? amp : one(T)
                data[end-i,j] *= right ? amp : one(T)
            end
        end
    end

    if form == :cosine
        # 1/2n without the π factor (since we call sinpi)
        omega_π = 1/T(2n)

        for j in axes(data, 2)
            @inbounds for i in 0:(n - 1)
                amp = sinpi(omega_π*T(i))
                data[begin+i,j] *= left ? amp : one(T)
                data[end-i,j] *= right ? amp : one(T)
            end
        end
    end

    nothing    
end
