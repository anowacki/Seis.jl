# Functions to filter traces
"""
    filt(f, t::AbstractTrace, args...) -> t′
    
Apply a filter `f` to a `Seis.Trace` `t` and return a modified copy.
"""
function DSP.filt(f, t::AbstractTrace, args...)
    t′ = deepcopy(t)
    DSP.filt!(trace(t′), f, trace(t), args...)
    t′
end

"""
    filt!(f, t::AbstractTrace, args...) -> t
    
Apply a filter `f` to a `Seis.Trace` `t` in-place and return the modified
original trace.
"""
function DSP.filt!(f, t::AbstractTrace, args...)
    # FIXME: This can be replaced with
    #     trace(t) .= DSP.filt(f, trace(t), args...)
    # because the RHS will be allocated before the LFS is filled in, and
    # there is therefore no aliasing.
    trace′ = deepcopy(trace(t))
    trace(t) .= DSP.filt(f, trace′, args...)
    t
end

"""
    filtfilt!(coef, t::AbstractTrace) -> t
    filtfilt(coef, t::AbstractTrace) -> t′
    
Apply a set of filter coefficients `coef` to the `Seis.Trace` `t` twice: both
forwards and backwards.
In the first form, update the trace in place and return it.
In the second form, return an updated copy.
"""
function filtfilt!(coef, t::AbstractTrace)
    trace′ = deepcopy(trace(t))
    trace(t) .= DSP.filtfilt(coef, trace′)
    t
end

"""
    bandpass!(t::Trace, f1, f2; poles=2, twopass=false) -> t
    bandpass(t::Trace, f1, f2; poles=2, twopass=false) -> t′
    
Apply a bandpass filter to the trace `t` between frequencies `f1` and `f2` in Hz.
In the first form, update the trace in place and return it.  In the second form,
return an updated copy.

Optionally specify the number of `poles` of the filter, and whether a two-pass
filter should be applied.  This doubles the effective number of poles, because
both a forward and reverse pass occur.  This has the advantage of preserving the
phase.

Specify the `kind` of filter by providing a kind from the `DSP.Filters` module.
If doing so, the `poles` keyword argument is not used and the number of poles,
ripple power, etc., should be specified when providing the filter kind.
"""
function bandpass!(t::AbstractTrace, f1, f2;
                   poles=2, twopass=false, kind=DSP.Butterworth(poles))
    filter = DSP.digitalfilter(DSP.Bandpass(f1, f2), kind; fs=1/t.delta)
    _apply_filter!(t, filter, twopass)
end
bandpass(t, args...; kwargs...) = bandpass!(deepcopy(t), args...; kwargs...)
@doc (@doc bandpass!) bandpass

"""
    bandstop!(t::Trace, f1, f2; poles=2, twopass=false, kind=DSP.Butterworth(poles)) -> t
    bandstop(t::Trace, f1, f2; poles=2, twopass=false, kind=DSP.Butterworth(poles)) -> t′

Apply a bandreject filter to the trace `t` with stop band between frequencies
`f1` and `f2` in Hz.
In the first form, update the trace in place and return it.  In the second form,
return an updated copy.

Optionally specify the number of `poles` of the filter, and whether a two-pass
filter should be applied.  This doubles the effective number of poles, because
both a forward and reverse pass occur.  This has the advantage of preserving the
phase.

Specify the `kind` of filter by providing a kind from the `DSP.Filters` module.
If doing so, the `poles` keyword argument is not used and the number of poles,
ripple power, etc., should be specified when providing the filter kind.
"""
function bandstop!(t::AbstractTrace, f1, f2;
                   poles=2, twopass=false, kind=DSP.Butterworth(poles))
    filter = DSP.digitalfilter(DSP.Bandstop(f1, f2), kind; fs=1/t.delta)
    _apply_filter!(t, filter, twopass)
end
bandstop(t, args...; kwargs...) = bandstop!(deepcopy(t), args...; kwargs...)
@doc (@doc bandstop!) bandstop

"""
    highpass!(t::Trace, f; poles=2, twopass=false, kind=DSP.Butterworth(poles)) -> t
    highpass(t::Trace, f; poles=2, twopass=false, kind=DSP.Butterworth(poles)) -> t′

Apply a highpass filter to the trace `t` with corner frequency `f` in Hz.
In the first form, update the trace in place and return it.  In the second form,
return an updated copy.

Optionally specify the number of `poles` of the filter, and whether a two-pass
filter should be applied.  This doubles the effective number of poles, because
both a forward and reverse pass occur.  This has the advantage of preserving the
phase.

Specify the `kind` of filter by providing a kind from the `DSP.Filters` module.
If doing so, the `poles` keyword argument is not used and the number of poles,
ripple power, etc., should be specified when providing the filter kind.
"""
function highpass!(t::AbstractTrace, f;
                   poles=2, twopass=false, kind=DSP.Butterworth(poles))
    filter = DSP.digitalfilter(DSP.Highpass(f), kind; fs=1/t.delta)
    _apply_filter!(t, filter, twopass)
end
highpass(t, args...; kwargs...) = highpass!(deepcopy(t), args...; kwargs...)
@doc (@doc highpass!) highpass

"""
    lowpass!(t::Trace, f; poles=2, twopass=false, kind=DSP.Butterworth(poles)) -> t
    lowpass(t::Trace, f; poles=2, twopass=false, kind=DSP.Butterworth(poles)) -> t′

Apply a lowpass filter to the trace `t` with corner frequency `f` in Hz.
In the first form, update the trace in place and return it.  In the second form,
return an updated copy.

Optionally specify the number of `poles` of the filter, and whether a two-pass
filter should be applied.  This doubles the effective number of poles, because
both a forward and reverse pass occur.  This has the advantage of preserving the
phase.

Specify the `kind` of filter by providing a kind from the `DSP.Filters` module.
If doing so, the `poles` keyword argument is not used and the number of poles,
ripple power, etc., should be specified when providing the filter kind.
"""
function lowpass!(t::AbstractTrace, f;
                  poles=2, twopass=false, kind=DSP.Butterworth(poles))
    filter = DSP.digitalfilter(DSP.Lowpass(f), kind; fs=1/t.delta)
    _apply_filter!(t, filter, twopass)
end
lowpass(t, args...; kwargs...) = lowpass!(deepcopy(t), args...; kwargs...)
@doc (@doc lowpass!) lowpass

"""
    _apply_filter!(t::Trace, filter, twopass) -> t

Apply a filter `filter` to the trace `t` and return the modified trace.

`filter` will usually be one of `Highpass`, `Lowpass`, `Bandpass` or `Bandstop`.

`kind` will usually be one of `Butterworth`, `Chebyshev1, `Chebyshev2`
or `Elliptic` from the `DSP.Filters` module.
"""
function _apply_filter!(t::AbstractTrace, filter, twopass)
    twopass ? filtfilt!(filter, t) : DSP.filt!(filter, t)
end
