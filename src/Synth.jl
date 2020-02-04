"""
# Seis.Synth

This module contains the following exported functions for creating synthetic `Trace`s:

- `gaussian1`: First-derivative Gaussian function.
- `heaviside`: Step function.
- `ricker`: Second-derivative Gaussian function.
- `sines`: Sum of an arbitrary number of sine waves.
- `spikes`: Single-sample spikes at arbitrary times in the trace.

*N.B.*  The main Seis module does not re-export these functions, so to bring them
into your namespace, do `using Seis.Synth`.
"""
module Synth

using ..Seis

export gaussian1, heaviside, ricker, sines, spikes

"""
    gaussian1(b, delta, n, ν, time=b+delta*n÷2; T=Float64, V=Vector{T}, P=Seis.Geographic{T}) -> t

Return a `Trace` `t` with start time `b` s, sampling interval `delta` s and
number of samples `n`, which has a synthetic signal consisting of a first-
derivative Gaussian centred at time `time` s with an approximate dominant
frequency `ν` Hz.

The maximum amplitude is 1.
"""
function gaussian1(b, delta, n, ν, time=b+delta*n÷2; T=Float64, V=Vector{T}, P=Seis.Geographic{T})
    ω = 2π*ν
    t = Trace{T,V,P}(b, delta, n)
    dtimes = times(t) .- time
    @. t.t = -ω * dtimes * exp(-ω^2*dtimes^2/2)
    t.t ./= maximum(abs, t.t)
    t
end

"""
    heaviside(b, delta, n, time=b+delta*n÷2; reverse=false, T=Float64, V=Vector{T}, P=Seis.Geographic{T}) -> t::Trace{T,V,S}

Return a `Trace` `t` with start time `b` s, sampling interval `delta` s and
number of samples `n`, which has a synthetic signal consisting of a step
from 0 before `time` s and 1 from `time` s onwards.

To reverse the function so that times on and before `time` are 1, use `reverse=true`.
"""
function heaviside(b, delta, n, time=b+delta*n÷2;
                   reverse=false, T=Float64, V=Vector{T}, P=Seis.Geographic{T})
    t = Trace{T,V,P}(b, delta, zeros(T, n))
    it = nearest_sample(t, time, inside=true)
    it === nothing && throw(ArgumentError("step time is outside the trace"))
    if reverse
        t.t[1:it] .= 1
    else
        t.t[it:end] .= 1
    end
    t
end

"""
    ricker(b, delta, n, ν, time=b+delta*n÷2; T=Float64, V=Vector{T}, P=Seis.Geographic{T}) -> t

Return a `Trace` `t` with start time `b` s, sampling interval `delta` s and
number of samples `n`, which has a synthetic signal consisting of a Ricker wavelet
centred at time `time` s with dominant frequency `ν` Hz.

The maximum amplitude is 1.
"""
function ricker(b, delta, n, ν, time=b+delta*n÷2; T=Float64, V=Vector{T}, P=Seis.Geographic{T})
    ω = 2π*ν
    t = Trace{T,V,P}(b, delta, n)
    dtimes = times(t) .- time
    @. t.t = (1 - ω^2*dtimes^2/2) * exp(-ω^2*dtimes^2/4)
    t.t ./= maximum(abs, t.t)
    t
end

"""
    sines(b, delta, n, νs, amps, ϕs=zeros(length(νs)); T=Float64, V=Vector{T}, P=Seis.Geographic{T}) -> t::Trace{T,V,S}

Return a trace object with start time `b`, sampling interval `delta` and number of
samples n, which has a synthetic signal consisting of `length(νs)` monochromatic
sine waves with frequencies `νs` in Hz, amplitudes `amps` and phases `ϕs` in radians.

No check is performed that any of the `νs` are below the Nyquist frequency determined
by `delta`.
"""
function sines(b, delta, n, νs, amps=ones(Float64, length(νs)),
               ϕs=zeros(Float64, length(νs)); T=Float64, V=Vector{T}, P=Seis.Geographic{T})
    length(νs) == length(amps) == length(ϕs) ||
        throw(ArgumentError("Arrays νs, amps and ϕs must be the same length"))
    t = Trace{T,V,P}(b, delta, zeros(T, n))
    for (ν, amp, ϕ) in zip(νs, amps, ϕs)
        t.t .+= amp*sin.(ν.*2π.*times(t) .- ϕ)
    end
    t
end

"""
    spikes(b, delta, n, times, amps=ones(length(times))) -> t::Trace

Return a trace `t` with single-sample spikes at `times` s, with amplitudes `amps`.
"""
function spikes(b, delta, n, times, amps=ones(length(times));
                T=Float64, V=Vector{T}, P=Seis.Geographic{T})
    length(times) == length(amps) ||
        throw(ArgumentError("arrays `times` and `amps` must be the same length"))
    t = Trace{T,V,P}(b, delta, zeros(T, n))
    for (time, amp) in zip(times, amps)
        it = nearest_sample(t, time, inside=true)
        it === nothing && error("spike time is outside the trace")
        t.t[it] = amp
    end
    t
end

end # module
