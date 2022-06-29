# Types and methods for Fourier-domain traces

"""
    AbstractFourierTrace <: AbstractRecording

Abstract type from which you should subtype if creating new types of
frequency-domain traces.  Note that this is a subtype of `AbstractData`,
meaning methods which work for `AbstractData` objects should work for
`AbstractFourierTrace` objects too.

# Interface

!!! note
    The formal interface for `AbstractFourierTrace`s is still a work in progress
    and may change in a minor version increment.

The following methods should be defined for all `AbstractFourierTrace`s `f`:

- `trace(f)`: Return the data for the Fourier trace.
- `frequencies(f)`: Return the frequency at each point of the Fourier series.
- `times(f)`: Return the time at each sample of the time domain trace
  corresponding to `f`.
- `starttime(f)`: The time of the first sample of the corresponding time
  domain trace.
- `nsamples(f)`: The number of samples in the time domain trace corresponding to `f`.
- `nfrequencies(f)`: The number of frequencies in the Fourier trace `f`.
- `Base.eltype(f)`: The element type of the frequency domain data samples.
"""
abstract type AbstractFourierTrace <: AbstractData end

"""
    FourierTrace <: AbstractFourierTrace

The Fourier transform of a [`Trace`](@ref).

`FourierTrace`s contain all the information contained in their corresponding
`Trace`, but with the addition of information on the original number of
samples in the time domain trace, and an element type which much be `Complex`.
Therefore, the fields `b`, `delta`, `evt`, `sta`, `picks`, and `meta` are all
part of the public API and can be accessed the same way as for `Trace`s.

The usual way to create `FourierTrace`s is by calling [`fft`](@ref) on a `Trace`.
To convert back to a `Trace`, call [`ifft`](@ref).

!!! note
    `FourierTrace`s share many fields with the `Trace` they came from, including
    the station and event, plus picks and metadata.  These fields are **not**
    copied across, and are instead references to the same `Station`, `Event`,
    and so on.  Therefore, any changes to the `FourierTrace` will be reflected
    in the corresponding `Trace`.  To avoid this, do `f = fft(deepcopy(t))`.

---
    FourierTrace{T,V,P}(b, delta, nsamples, data, evt, sta, picks, meta)

Create a new `FourierTrace` which corresponds to a `Trace` with start time
`b` s, sampling interval `delta` s and `nsamples` data points.  `data`
contains a set of Fourier coefficients, starting at 0 frequency up to
the Nyquist frequency.  (The coefficient normalisation is defined to be
the same as that used by FFTW.)  `evt`, `sta`, `picks` and `meta` should
be taken from the original `Trace`.

!!! note
    One usually creates `FourierTrace`s by calling `fft` on a `Trace`.
    Users are advised to use the keyword constructor (below) if
    creating new `FourierTrace`s from scratch.

---
    FourierTrace(; b, delta, data, nsamples=(2*length(data) - 1), evt=Event(), sta=Station(), picks=nothing, meta=nothing)

Create a `FourierTrace` directly from a set of Fourier coefficients.
Users will usually construct a `FourierTrace` by calling [`fft`](@ref)
on a `Trace` instead of using this constructor.

# Keyword arguments
- `b`: The starting time in s of the original recording this frequency
  domain trace represents.
- `delta`: The original sampling interval in s of the equivalent time
  domain trace.
- `data`: An `AbstractArray{<:Complex}` containing the set of Fourier
  coefficients for this frequency domain trace.  Note that this is a
  'one-sided' set, where the first index corresponds to 0 Hz and the final
  index is the Nyquist frequency.  This is because `Trace`s represent real
  quantities.
- `nsamples`: The number of time domain samples in the origin time
  domain trace which this frequency domain trace represents.  This allows
  the `FourierTrace` to be converted back to the original `Trace` with no
  loss of the number of points.
- `evt`: An `Event` which defines the source and origin time for the data.
  Normally this should be taken from the `Trace` being used to construct
  the `FourierTrace`.
- `sta`: A `Station` which defines the recording station for the data.
  Normally this should be taken from the `Trace` being used to construct
  the `FourierTrace`.

---
See also: [`Trace`](@ref), [`fft`](@ref), [`ifft`](@ref).
"""
mutable struct FourierTrace{T<:AbstractFloat,
                            V<:AbstractVector{<:Complex{<:AbstractFloat}},
                            P<:Position{T}
                            } <: AbstractFourierTrace
    b::T
    # Frequency spacing in Hz, not  sampling interval in s
    delta::T
    nsamples::Int64
    data::V
    evt::Event{T,P}
    sta::Station{T,P}
    picks::SeisDict{Union{Symbol,Int},Pick{T}}
    meta::SeisDict{Symbol,Any}

    function FourierTrace{T,V,P}(b, delta, nsamples, data, evt, sta, picks, meta
            ) where {T,V,P}
        delta > 0 || throw(ArgumentError("delta cannot <= 0"))
        nsamples >= 0 || throw(ArgumentError("nsamples cannot be negative"))
        new{T,V,P}(b, delta, nsamples, data, evt, sta, picks, meta)
    end
end

"""
    FourierTrace(; b, delta, data, nsamples=(2*length(data) - 1), evt=Event(), sta=Station(), picks=nothing, meta=nothing)

Create a `FourierTrace` directly from a set of Fourier coefficients.
Users will usually construct a `FourierTrace` by calling [`fft`](@ref)
on a `Trace` instead of using this constructor.

# Keyword arguments
- `b`: The starting time in s of the original recording this frequency
  domain trace represents.
- `delta`: The original sampling interval in s of the equivalent time
  domain trace.
- `data`: An `AbstractArray{<:Complex}` containing the set of Fourier
  coefficients for this frequency domain trace.  Note that this is a
  'one-sided' set, where the first index corresponds to 0 Hz and the final
  index is the Nyquist frequency.  This is because `Trace`s represent real
  quantities.
- `nsamples`: The number of time domain samples in the origin time
  domain trace which this frequency domain trace represents.  This allows
  the `FourierTrace` to be converted back to the original `Trace` with no
  loss of the number of points.
- `evt`: An `Event` which defines the source and origin time for the data.
  Normally this should be taken from the `Trace` being used to construct
  the `FourierTrace`.
- `sta`: A `Station` which defines the recording station for the data.
  Normally this should be taken from the `Trace` being used to construct
  the `FourierTrace`.
"""
function FourierTrace{T,V,P}(;
        b, delta, data, nsamples=(2*length(data) - 1),
        evt=Event{T,P}(), sta=Station{T,P}(),
        picks=SeisDict{Union{Symbol,Int},Pick{T}}(),
        meta=SeisDict{Symbol,Any}()
    ) where {T,V,P}
    FourierTrace{T,V,P}(b, delta, nsamples, data, evt, sta, picks, meta)
end

FourierTrace{T,V}(; kwargs...) where {T,V} = FourierTrace{T,V,Geographic{T}}(; kwargs...)
FourierTrace{T}(; kwargs...) where {T} =
    FourierTrace{T,Vector{Complex{T}},Geographic{T}}(; kwargs...)
FourierTrace(; kwargs...) =
    FourierTrace{Float64,Vector{Complex{Float64}},Geographic{Float64}}(; kwargs...)

Base.eltype(f::AbstractFourierTrace) = eltype(trace(f))


@eval begin
    function Base.hash(f::FourierTrace, h::UInt)
        $([:(h = hash(f.$f, h)) for f in fieldnames(FourierTrace)]...)
        h
    end

    function Base.:(==)(a::FourierTrace, b::FourierTrace)
        all(isequal(getfield(a, f), getfield(b, f)) for f in $(fieldnames(FourierTrace)))
    end
end

# Get an array of values from an array of Fourier traces
function Base.getproperty(t::AbstractArray{<:AbstractFourierTrace}, f::Symbol)
    if f === :b || f === :delta || f === :nsamples || f === :data ||
            f === :evt || f === :sta || f === :picks || f === :meta
        getfield.(t, f)
    else
        getfield(t, f)
    end
end

# Set an array of traces with a single value
function Base.setproperty!(t::AbstractArray{<:AbstractFourierTrace}, f::Symbol, val)
    if f === :b || f === :delta || f === :meta || f === :evt || f == :sta || f === :picks
        for tt in t
            setproperty!(tt, f, val)
        end
    else
        # Fallback in case of desired dot-access to fields of an array type
        setfield!(t, f, val)
    end
end

# Set an array of traces with an array of values
function Base.setproperty!(t::AbstractArray{<:AbstractFourierTrace}, f::Symbol, val::AbstractArray)
    length(t) == length(val) || throw(DimensionMismatch())
    for (tt, vv) in zip(t, val)
        setproperty!(tt, f, vv)
    end
end

Base.broadcastable(f::FourierTrace) = Ref(f)

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
