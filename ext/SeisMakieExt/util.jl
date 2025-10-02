"""
    _decimation_value(trace, shift, t1, t2, max_samples) -> n

Return the decimation value `n` required to reduce the number
of samples in `trace` between `t1` and `t2`, with a time shift
of `shift` s, to max_samples
"""
function _decimation_value(trace::Seis.AbstractTrace, shift, t1, t2, max_samples)
    npts = Seis.nsamples(trace, t1 + shift, t2 + shift)
    max(1, round(Int, npts/max_samples))
end
function _decimation_value(times, t1, t2, max_samples)
    npts = sum(x -> t1 <= x <= t2, times)
    max(1, round(Int, npts/max_samples))
end

"Return the maximum absolute value over all traces"
_max_abs_value(ts...) = maximum(t -> maximum(abs, Seis.trace(t)), ts)

function _kwargs_deprecation(old_name, new_name, kwargs_old, kwargs_new)
    if !isnothing(kwargs_old)
        @warn("$old_name is deprecated in favour of $new_name and will be removed in a future version=")
        kwargs_old
    else
        kwargs_new
    end
end
