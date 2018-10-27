"""
    endtime(t) -> time

Return the end `time` of trace `t` in seconds.
"""
endtime(t::AbstractTrace) = t.b + (nsamples(t) - 1)*t.delta

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
