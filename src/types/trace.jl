"""
    AbstractTrace <: AbstractData

Abstract type from which you should subtype if creating new types of traces.
`AbstractTrace`s are time-domain recordings (or synthetics) at a single
channel.

# Interface

!!! note
    The formal interface for `AbstractTrace`s is still a work in progress and may change with
    a minor version increment.

The following methods should be defined for all `AbstractTrace`s `t`:

- `trace(t)`: Return the data for the trace.
- `times(t)`: Return the time at each sample of `t`.
- `starttime(t)`: The time of the first sample.
- `nsamples(t)`: The number of samples in `t`.
- `Base.eltype(t)`: The element type of the data samples.
- `t.evt`: Return the `Event` associated with this trace.
- `t.sta`: Return the `Station` at which this trace was recorded.
- `t.meta`: Return a `SeisDict{Symbol,Any}` into which metadata may be placed.
"""
abstract type AbstractTrace end

"""
    Trace

Evenly-sampled time series recorded at a single seismic station.  The start time
of the trace, in s, is in the `b` property, whilst the sampling interval, in s,
is `delta`.
The trace itself is accessed using the `trace` method, like `trace(t)`.

All `Trace`s are relative to the event time `evt.time` if it is defined, regardless
of what the event is.  For example, `evt.time` could be the origin time of an
earthquake, or a picked arrival time.

The trace then contains information about an associated `Event` in `evt` and the
`Station` in `sta`.  `picks` holds a dictionary which contains pairs of pick times
relative to the origin and names (which can be `missing`).  Access picks with the
`picks` method, and add picks with `add_pick!`.

The `meta` `Dict` holds any other information about the trace.

If the event `time` is set, then the trace beginning time `b` is relative to this.

Find the trace start time relative to the origin time using [`starttime`](@ref).
The absolute start time and date, if an origin time is set, is given by
[`startdate`](@ref).
"""
mutable struct Trace{T<:AbstractFloat,V<:AbstractVector{<:AbstractFloat},P<:Position{T}} <: AbstractTrace
    b::T
    delta::T
    t::V
    evt::Event{T,P}
    sta::Station{T,P}
    picks::SeisDict{Union{Symbol,Int},Pick{T}}
    meta::SeisDict{Symbol,Any}
    function Trace{T,V,P}(b, delta, t, evt, sta, picks, meta) where {T,V,P}
        delta > 0 || throw(ArgumentError("delta cannot be <= 0"))
        new(b, delta, t, evt, sta, picks, meta)
    end
end

"""
    Trace(b, delta, t::AbstractVector) -> trace::Trace{$DEFAULT_FLOAT,Vector{$DEFAULT_FLOAT},Seis.Geographic{$DEFAULT_FLOAT}}

Create a `Trace` with starting time `b` s, sampling interval `delta` s and
an `AbstractVector` `t` of values for the trace.  The default precision for
the type is `$DEFAULT_FLOAT`.  By default, the trace's event and station
are in geographic coordinates.

    Trace(b, delta, n::Integer) -> trace::Trace{$DEFAULT_FLOAT,Vector{$DEFAULT_FLOAT},Seis.Geographic{$DEFAULT_FLOAT}}

Create a new `Trace` with uninitialised data of length `n` samples.

    Trace{T,V,P}(args...) -> trace::Trace{T,V,P}
    Trace{T,V}(args...) -> trace::Trace{T,V,Seis.Geographic{T}}
    Trace{T}(args...) -> trace::Trace{T,Vector{T},Seis.Geographic{T}}

Construct traces with non-default number and data types.  In the third
form, the data type defaults to `Vector{$DEFAULT_FLOAT}`, whilst in
both the second and third forms, `P` defaults to `Seis.Geographic{T}`.
"""
Trace{T,V,P}(b, delta, t::AbstractVector) where {T,V,P} =
    Trace{T,V,P}(b, delta, t, Event{T,P}(), Station{T,P}(), Dict(), Dict())
Trace{T,V,P}(b, delta, n::Integer) where {T,V,P} = Trace{T,V,P}(b, delta, Vector{T}(undef, n))
Trace{T,V}(args...) where {T,V} = Trace{T,V,Geographic{T}}(args...)
Trace{T}(args...) where T = Trace{T,Vector{T},Geographic{T}}(args...)
Trace(b, delta, t_or_n) = Trace{DEFAULT_FLOAT,Vector{DEFAULT_FLOAT},Geographic{DEFAULT_FLOAT}}(b, delta, t_or_n)

"""
    CartTrace

Alias for `Trace` where `Event` and `Station` coordinates are
`Seis.Cartesian` rather than `Seis.Geographic`
"""
const CartTrace{T,V} = Trace{T, V, Cartesian{T}}

"""
    CartTrace

Create a new `Trace` whose `Event` and `Station` are in Cartesian coordinates.

See [`Trace`](@ref) for more details of the different construction methods.

---
    CartTrace{T,V}(b, delta, t::AbstractVector) -> trace
    CartTrace{T}(args...)

Construct a `Trace` with Cartesian coordinates for the `Event` and `Station`
with non-default number and data types.  See [`Trace{T,V}`](@ref)
for details.
"""
CartTrace{T}(args...) where T = Trace{T, Vector{T}, Cartesian{T}}(args...)
CartTrace(args...) = Trace{DEFAULT_FLOAT, Vector{DEFAULT_FLOAT},
          Cartesian{DEFAULT_FLOAT}}(args...)

const TRACE_FIELDS = fieldnames(Trace)

# Get an array of values from an array of traces
function Base.getproperty(t::AbstractArray{<:Trace}, f::Symbol)
    if f === :b || f === :delta || f === :meta || f === :evt || f === :sta || f === :picks
        getfield.(t, f)
    else
        getfield(t, f)
    end
end

# Set an array of traces with a single value
function Base.setproperty!(t::AbstractArray{<:Trace}, f::Symbol, val)
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
function Base.setproperty!(t::AbstractArray{<:Trace}, f::Symbol, val::AbstractArray)
    length(t) == length(val) || throw(DimensionMismatch())
    for (tt, vv) in zip(t, val)
        setproperty!(tt, f, vv)
    end
end

Base.propertynames(t::AbstractArray{<:Trace}, private=false) = fieldnames(eltype(t))

Base.:(==)(t1::Trace, t2::Trace) =
    all(x -> isequal(x[1], x[2]), (getfield.((t1, t2), f) for f in TRACE_FIELDS))


# Treat single objects as scalars in broadcasting
Base.broadcastable(t::Union{Trace,Event,Station,Position,Pick}) = Ref(t)

# Element type of trace
Base.eltype(::Trace{T,V,P}) where {T,V,P} = eltype(V)
