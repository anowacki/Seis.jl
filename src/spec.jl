# Methods for Fourier-domain traces

"""
    fft(t::Trace) -> f::FourierTrace

Convert the trace `t` into its equivalent frequency domain trace `f` by
performing a Fourier transform.  The object returned is a [`FourierTrace`](@ref).

# Example
```
julia> t = Trace(0, 0.01, 1000); trace(t) .= sin.(2π.*times(t)); # Sine wave of frequency 1 Hz

julia> f = fft(t);

julia> indmax = argmax(abs.(trace(f))); # Get index of maximum power

julia> frequencies(f)[indmax] # Maximum frequency in Hz is 1 as expected
1.0
```

See also: [`ifft`](@ref).
"""
function fft(t::Trace{T,<:AbstractVector{TV},P}) where {T,TV,P}
    data = FFTW.rfft(trace(t))
    delta = 1/(nsamples(t)*t.delta)
    FourierTrace{T,Vector{Complex{TV}},P}(starttime(t), delta, nsamples(t),
        data, t.evt, t.sta, t.picks, t.meta)
end

"""
    ifft(t::FourierTrace[, d=2*nfrequencies(f) - 2]) -> t::Trace

Convert the frequency domain [`FourierTrace`](@ref) `f` back into its
equivalent time domain [`Trace`](@ref) `t`.

Because `FourierTrace`s are usually constructed by calling [`fft`](@ref)
on a `Trace`, the original number of time samples is kept in `f` and
`t` therefore contains the same number of samples as before so long as
the raw data in `f` has not been shortened or lengthened.  If it has,
then it is assumed that `t` should have an even number of points.  If
instead `t` should have an odd number of point, pass `d=npts`, where `npts`
is the number of points needed.  `d` must be either one or two less than
double the number of frequency points in `t`.

# Example

Showing that taking the inverse Fourier transform of the Fourier domain trace
`f` gets us back to `t`:
```
julia> t = sample_data();

julia> f = fft(t);

julia> t′ = ifft(f);

julia> isapprox(trace(t), trace(t′))
true
```

A crude method to upsample a trace:
```
julia> t = Trace(0, 1, [0, 1, 0, -1, 0]);

julia> f = fft(t);

julia> append!(trace(f), zeros(nfrequencies(f)));

julia> trace(ifft(f))
10-element Vector{Float64}:
  0.0
  0.447213595499958
  1.0
  0.894427190999916
 -1.7763568394002506e-16
 -0.894427190999916
 -1.0
 -0.44721359549995787
  1.7763568394002506e-16
  0.0
```

See also: [`fft`](@ref).
"""
function ifft(f::FourierTrace{T,<:AbstractVector{<:Complex{TV}},P},
        d::Union{Integer,Nothing}=nothing) where {T,TV,P}
    # Use original number of points if the length of the frequency coefficients
    # array has not been changed
    n = if f.nsamples in 2*nfrequencies(f) .- (1, 2)
        f.nsamples
    # Otherwise, we can't know how long the time-domain trace is meant
    # to be and so assume we want an even number of points unless we are told
    # how many we should have
    elseif d === nothing
        2*nfrequencies(f) - 2
    else
        if d == 2*nfrequencies(f) - 1 || d == 2*nfrequencies(f) - 2
            d
        else
            throw(ArgumentError("d must be either 2nf - 2 or 2nf - 1, " *
                "where nf is the number of frequency points"))
        end
    end

    # Scale the amplitude appropriately in case we have changed the
    # spectrum length
    scale = n/f.nsamples
    # Ensure we pass `d` as an `Int` (not `Int64`) on x86, otherwise there
    # is no method to construct a `FFTW.FakeArray`.
    data = FFTW.irfft(trace(f), Int(n)).*scale

    delta = 1/(2*(nfrequencies(f) - 1)*f.delta)
    Trace{T,Vector{TV},P}(f.b, delta, data, f.evt, f.sta, f.picks, f.meta)
end

"""
    nsamples(f::AbstraceFourierTrace[; even::Bool])

For a `FourierTrace` `f`, return the number of samples in the equivalent
time domain trace.

The Fourier trace `f` records the number of points in the original trace
used to create it with a call to [`fft`](@ref), however if the length of
the trace has been changed, the number of samples of the time-domain
equivalent of `f` may be odd or even.  In this case, pass either `even=true`
to return the even number of sample, or `even=false` for the odd.

# Example
```
julia> t = sample_data(); nsamples(t)
1000

julia> f = fft(t);

julia> nsamples(f)
1000
```

See also: [`nfrequencies`](@ref)
"""
function nsamples(f::AbstractFourierTrace; even::Union{Nothing,Bool}=nothing)
    nfreq = nfrequencies(f)
    if f.nsamples in (2nfreq - 1, 2nfreq - 2)
        if even === nothing
            return f.nsamples
        else
            return even ? 2nfreq - 2 : 2nfreq - 1
        end
    else
        return (even === nothing || even) ? 2nfreq - 2 : 2nfreq - 1
    end
end

"""
    frequencies(f)

Return the frequency in Hz of each data point in the frequency domain
trace `f`.

# Example
```
julia> f = fft(sample_data());

julia> frequencies(f)
0.0f0:0.1f0:50.0f0
```
"""
frequencies(f::AbstractFourierTrace) = f.delta.*(0:(nfrequencies(f) - 1))

"""
    nfrequencies(f)

Return the number of frequency points in the frequency domain trace `f`.

# Example
```
julia> f = fft(sample_data());

julia> nfrequencies(f)
501
```

See also: [`nsamples`](@ref)
"""
nfrequencies(f) = length(trace(f))

"""
    trace(f::FourierTrace)

Return the data for a [`FourierTrace`](@ref) object.
"""
trace(f::AbstractFourierTrace) = f.data

"""
    starttime(f::FourierTrace)

Return the time of the first sample of the equivalent time-domain recording
which would be obtained by [`ifft`](@ref ifft(::FourierTrace))
"""
starttime(f::AbstractFourierTrace) = f.b
