"""
    angle_difference(α, β, degrees::Bool=true) -> Δ

Return the angular difference `Δ` between two angles `α` and `β`, in the direction
α → β (i.e., Δ = (β - α)).  This means that Δ is positive if the direction β is
clockwise of α, and negative otherwise.

Angles are assumed to be degreees unless `degrees` is `false`.
"""
function angle_difference(a, b, degrees::Bool=true)
    T = float(promote_type(typeof(a), typeof(b)))
    whole_circ = degrees ? T(360) : T(2)*π
    half_circ = whole_circ/2
    mod(b - a + half_circ, whole_circ) - half_circ
end

"""
    _angle_tol(x) -> tol

Return an appropriate tolerance for `x` (a floating point type, `AbstractTrace` or
`Station`) to use in comparisons of angular quantities in °.  Angles which are `tol`°
apart can be considered identical.
"""
_angle_tol(T::DataType) = √eps(T)
_angle_tol(::Station{T}) where T = _angle_tol(T)
_angle_tol(t::AbstractTrace) = _angle_tol(t.sta)

"""
    dates(t) -> date_range

Return a `date_range` which contains the dates for each sample of `t`, so long
as `t.evt.time` is defined.  If not, an error is thrown.

N.B.  This function assumes that the sampling interval `t.delta` is representable
as an integer number of milliseconds, and rounds it accordingly.  `Dates.DateTime`s
have precision of 1 ms.  An error is thrown if `t.delta < 1e-3` s.

# Example
```
julia> t = sample_data();

julia> dates(t)
1981-03-29T10:39:06.66:10 milliseconds:1981-03-29T10:39:16.65
```

See also: [`times`](@ref).
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

# Example
```
julia> t = sample_data(); t.evt.time
1981-03-29T10:38:14

julia> startdate(t)
1981-03-29T10:39:06.66
```
"""
startdate(t::AbstractTrace) = ((b, delta) = _check_date_b_delta(t); t.evt.time + b)

"""
    enddate(t) -> date

Return the `date` of the last sample of the trace `t`.

N.B.  This function assumes that the sampling interval `t.delta` is representable
as an integer number of milliseconds, and rounds it accordingly.  `Dates.DateTime`s
have precision of 1 ms.  An error is thrown if `t.delta < 1e-3` s.

# Example
```
julia> t = sample_data(); t.evt.time
1981-03-29T10:38:14

julia> enddate(t)
1981-03-29T10:39:16.65
```
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

# Example
```
julia> t = Trace(-3, 0.01, rand(20)) # Set start time to -3;

julia> starttime(t)
-3.0
```
"""
starttime(t::AbstractTrace) = t.b

"""
    endtime(t) -> time

Return the end `time` of trace `t` in seconds.

# Example
```
julia> t = Trace(5, 1, 3); # 3 samples at 1 Hz, starting at 5 s

julia> endtime(t)
7.0
```
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

# Example
```
julia> using Dates

julia> t = sample_data(); t.picks
Seis.SeisDict{Union{Int64, Symbol},Seis.Pick{Float32}} with 2 entries:
  :F => Seis.Pick{Float32}(time=60.980003, name=missing)
  :A => Seis.Pick{Float32}(time=53.670002, name=missing)

julia> t.evt.time
1981-03-29T10:38:14

julia> origin_time!(t, t.evt.time + Second(1)); t.evt.time
1981-03-29T10:38:15

julia> t.picks
Seis.SeisDict{Union{Int64, Symbol},Seis.Pick{Float32}} with 2 entries:
  :F => Seis.Pick{Float32}(time=59.980003, name=missing)
  :A => Seis.Pick{Float32}(time=52.670002, name=missing)
```
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
            t.picks[key] = Seis.Pick{eltype(t)}(time-Δb, name)
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
    _has_azimuth(s::Station{T}, azi; tol=_angle_tol(T)) where T -> ::Bool

Return `true` if the azimuth of station `s` is within `tol` of `azi`.
Angles `azi` and `tol` are in degrees.
"""
function _has_azimuth(s::Station{T}, azi; tol=_angle_tol(T)) where T
    isapprox(abs(angle_difference(s.azi, azi)), zero(T), atol=tol)
end

"""
    is_east(s::Station; tol) -> ::Bool
    is_east(t::AbstractTrace; tol) -> ::Bool

Return `true` if the trace `t` or station `s` is horizontal and points to
the east.

The azimuth and inclination of the trace is compared to east and the horizontal
within a tolerance of `tol`°.  The default is set to be appropriate for the
floating-point type used for the station or trace, but can be overridden by
passing a comparison to `tol`.

See also: [`is_north`](@ref), [`is_vertical`](@ref)
"""
function is_east(s::Station{T}; tol=_angle_tol(T)) where T
    is_horizontal(s, tol=tol) && _has_azimuth(s, 90, tol=tol)
end
is_east(t::Trace{T}; tol=_angle_tol(T)) where T = is_east(t.sta, tol=tol)

"""
    is_north(s::Station{T}; tol) where T -> ::Bool
    is_north(t::AbstractTrace; tol) -> ::Bool

Return `true` if the trace `t` is horizontal and points to the north.

The azimuth and inclination of the trace is compared to east and the horizontal
within a tolerance of `tol`°.  The default is set to be appropriate for the
floating-point type used for the station or trace, but can be overridden by
passing a comparison to `tol`.

See also: [`is_east`](@ref), [`is_vertical`](@ref)
"""
function is_north(s::Station{T}; tol=_angle_tol(T)) where T
    is_horizontal(s, tol=tol) && _has_azimuth(s, 0, tol=tol)
end
is_north(t::Trace{T}; tol=_angle_tol(T)) where T = is_north(t.sta, tol=tol)

"""
    is_horizontal(s::Station; tol)
    is_horizontal(t::AbstractTrace; tol) -> ::Bool

Return `true` if the trace `t` is horizontal (i.e., its inclination
is 90° from the vertical), and `false` otherwise.

The inclination of the trace is compared to the horizontal
within a tolerance of `tol`°.  The default is set to be appropriate for the
floating-point type used for the station or trace, but can be overridden by
passing a comparison to `tol`.

# Examples
```
julia> s = Station(azi=0, inc=90);

julia> is_horizontal(s)
true

julia> t = sample_data();

julia> is_horizontal(t)
false

julia> t.sta.inc
0.0f0
```

See also: [`is_vertical`](@ref), [`is_east`](@ref), [`is_north`](@ref).
"""
is_horizontal(s::Station{T}; tol=_angle_tol(T))  where T = isapprox(s.inc, 90, atol=tol)
is_horizontal(t::Trace{T}; tol=_angle_tol(T)) where T = is_horizontal(t.sta, tol=tol)

"""
    is_vertical(s::Station{T}; tol=eps(T)) where T
    is_vertical(t::AbstractTrace; tol=eps(eltype(trace(t)))) -> ::Bool

Return `true` if the trace `t` is vertical (i.e., its inclination
is 0°), and `false` otherwise.

The inclination of the trace is compared to the vertical
within a tolerance of `tol`°.  The default is set to be appropriate for the
floating-point type used for the station or trace, but can be overridden by
passing a comparison to `tol`.

# Examples
```
julia> s = Station(azi=0, inc=90);

julia> is_vertical(s)
false

julia> t = sample_data();

julia> is_vertical(t)
true
```

See also: [`is_horizontal`](@ref).
"""
is_vertical(s::Station{T}; tol=_angle_tol(T))  where T = isapprox(s.inc, 0, atol=tol)
is_vertical(t::Trace{T}; tol=_angle_tol(T)) where T = is_vertical(t.sta, tol=tol)

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

# Examples
```
julia> t = Trace(0, 1, rand(5)); # Trace starting at 0 s, 1 Hz sampling

julia> nearest_sample(t, 2)
3

julia> nearest_sample(t, -1)

julia> nearest_sample(t, -1, inside=false)
1
```
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

# Example
```
julia> using Dates: DateTime, Second

julia> t = sample_data();

julia> nearest_sample(t, DateTime(1981, 03, 29, 10, 39, 7))
35

julia> nearest_sample(t, startdate(t) - Second(10)) # 10 s before the first sample

julia> nearest_sample(t, startdate(t) -  Second(10), inside=false)
1
```
"""
function nearest_sample(t::AbstractTrace, datetime::DateTime; inside=true)
    ismissing(t.evt.time) && error("trace does not have origin time set")
    d = dates(t)
    if inside
        d[1] <= datetime <= d[end] || return nothing
    end
    datetime <= d[1] && return 1
    datetime >= d[end] && return nsamples(t)
    argmin(abs.(d .- datetime))
end

"""
    nsamples(t) -> n

Return the number of samples `n` in a trace `t`.

# Example
```
julia> data = rand(4);

julia> t = Trace(0, 1, data);

julia> nsamples(t)
4
```
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
    ib = b == -Inf ? 1 : clamp(ceil(Int, (b - tb)/delta) + 1, 1, n)
    ie = e == Inf ? n : clamp(floor(Int, (e - tb)/delta) + 1, 1, n)
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

# Example
```
julia> t = sample_data();

julia> times(t)
52.66f0:0.01f0:62.65f0
```

See also: [`dates`](@ref).
"""
times(t::AbstractTrace) = t.b .+ (0:(nsamples(t) - 1)).*t.delta

"""
    trace(t) -> data

Return an array `data` containing the values of the `Trace` `t` at each
sampling point.  `data` is now a variable bound to `t`s values, and changing
`data` will change `t`.  `trace(t)` may itself also be modified and the trace
will be updated.

!!! note
    The value returned by `trace` is a variable bound to an internal field
    of the trace.  Therefore, assigning another value to `trace(t)` or
    `data` will **not** update the values in `t`.  Instead, update the values
    in-place using the `.` operator (like `data .= 1`).  See examples below.

!!! note
    The underlying data array holding the trace can be rebound by assigning to
    the trace's field `t`, but this is unsupported and may break in future.

# Examples

Retrieving the data values for a trace, and modifying the first value.
```
julia> t = sample_data();

julia> data = trace(t)
1000-element Array{Float32,1}:
 -0.09728001
 -0.09728001
  ⋮
 -0.0768
 -0.0768

julia> data[1] = 0;

julia> trace(t)
1000-element Array{Float32,1}:
  0.0
 -0.09728001
  ⋮
 -0.0768
 -0.0768
```

Setting the data values for a new synthetic trace.
```
julia> t = Trace(0, 0.01, 1000); # 1000-point, 100 Hz trace with random data

julia> trace(t) .= sin.(π.*times(t));

julia> trace(t)
1000-element Array{Float64,1}:
  0.0
  0.03141075907812829
  ⋮
 -0.06279051952931425
 -0.031410759078131116
```
"""
trace(t::AbstractTrace) = t.t

"""
    traces_are_orthogonal(t1::Trace, t2::Trace; tol=eps()) -> ::Bool

Return `true` if the two traces `t1` and `t2` have component azimuths
90° apart.  Set the tolerance of the comparison with `tol`.
"""
traces_are_orthogonal(t1::AbstractTrace, t2::AbstractTrace;
                      tol=max(_angle_tol(t1), _angle_tol(t1))) =
    isapprox(abs(angle_difference(t1.sta.azi, t2.sta.azi)), 90, atol=tol)
