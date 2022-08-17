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
