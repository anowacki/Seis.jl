"""
    Response

Abstract supertype of all response types, which describe the
transfer function between different stages in the recording of
a signal using a seismometer.  For example, a response may
represent the phase and amplitude shift between the ground
motion and the raw electrical signal internal to an instrument,
or how an antialiasing filter is then applied to that signal
before being digitised.
"""
abstract type Response end

"""
    PolesZerosResponse{T} <: Response

Type of instrument response which can be represented by a set
of complex poles and zeros, plus a real gain.

The transfer function ``H(\\omega)`` is given by

```math
H(\\omega) = G \\prod_{j=1}^M (i\\omega - z_j) / \\prod_{k=1}^N (i\\omega - p_k)``
```

where ``\\omega`` is the angular frequency in radians per second,
``z_j`` is the ``j``th zero, ``p_k`` is the ``k``th pole, and there
are ``M`` zeros and ``N`` poles.  ``G`` is the scalar gain factor
and ``i`` is the imaginary unit.

!!! note
    StationXML and SEED allow for poles and zeros to be in either
    Hz or radians/second.  `PolesZerosResponse`s are always in
    radians/second in Seis.jl.

# Accessible fields
The following fields are part of the public API for `PolesZerosResponse`
and should be accessed directly if need be:
- `poles::Vector{Complex{T}}`: Vector of complex poles.  May be empty.
- `zeros::Vector{Complex{T}}`: Vector of complex zeros.  May be empty.
- `gain::T`: Real gain of system.

---

    PolesZerosResponse(poles, zeros, gain)
    PolesZerosResponse(; poles, zeros, gain)

Construct a `PolesZerosResponse` from a set of `poles`, `zeros`
and a `gain`.  In the second form, keyword arguments may be used but
all must be present.

---

    PolesZerosResponse{T}(args...) -> ::PolesZerosResponse{T}

Construct a `PolesZerosResponse` with element type `T`.  In the other
forms, the element type is inferred from the arguments
"""
struct PolesZerosResponse{T<:Real} <: Response
    poles::Vector{Complex{T}}
    zeros::Vector{Complex{T}}
    gain::T
end

function PolesZerosResponse(poles::AbstractVector{<:Complex},
        zeros::AbstractVector{<:Complex}, gain)
    T = float(promote_type(real(eltype(poles)), real(eltype(zeros)), typeof(gain)))
    PolesZerosResponse{T}(poles, zeros, gain)
end    
PolesZerosResponse(; poles, zeros, gain) = PolesZerosResponse(poles, zeros, gain)

"""
    PolesZerosResponse(poles, nzeros::Integer, gain)

Constract a `PolesZerosResponse` when there are `nzeros` zeros
which are on the origin of the complex plane; i.e., they are
actually 0.
"""
function PolesZerosResponse(poles, nzeros::Integer, gain)
    zeros = zeros(eltype(poles), nzeros)
    PolesZerosResponse(poles, zeros, gain)
end

Base.broadcastable(paz::PolesZerosResponse) = Ref(paz)

function Base.show(io::IO, ::MIME"text/plain", paz::PolesZerosResponse{T}) where T
    println(io, "PolesZerosResponse{$T}:")
    println(io, "  poles: $(length(paz.poles)): $(paz.poles)")
    println(io, "  zeros: $(length(paz.zeros)): $(paz.zeros)")
    print(io, "   gain: $(paz.gain)")
end

"""
    frequency_response(response::Response, frequency) -> resp::Complex

Compute the transfer function of the `response` at `frequency` Hz.
The returned value `resp` is complex; to find the response amplitude,
do `abs(resp)`, and to find the phase response do `angle(resp)`.

# Example
Find the amplitude and phase (radians) response of a SmartSolo IGU16-HR
geophone at 100 Hz:
```
julia> paz = PolesZerosResponse([-22.211059 + 22.217768im, -22.211059 - 22.217768im], 2, 76.7)
PolesZerosResponse{Float64}:
  poles: 2: ComplexF64[-22.211059 + 22.217768im, -22.211059 - 22.217768im]
  zeros: 2: ComplexF64[0.0 + 0.0im, 0.0 + 0.0im]
   gain: 76.7

julia> resp = frequency_response(paz, 100)
-2.940646252079005 + 0.8662654098804758im

julia> amp, phase = abs(resp), angle(resp)
(76.69981822381358, 0.07075886038970215)
```
"""
function frequency_response(paz::PolesZerosResponse{T}, frequency) where T
    ω = T(2*frequency)*π
    iω = im*ω
    response = paz.gain*prod(z -> iω - z, paz.zeros; init=one(iω)) /
                        prod(p -> iω - p, paz.poles; init=one(iω))
end

"Union containing all response types"
const ResponseTypes{T} = Union{PolesZerosResponse{T}}

# Set of internal Symbols representing different response units
const UNITS_M = :m
const UNITS_M_S = :m_s
const UNITS_M_S2 = :m_s2
const UNITS_V = :V
const UNITS_COUNTS = :counts


struct ResponseStage{T}
    response::ResponseTypes{T}
    input_units::Symbol
    output_units::Symbol
    gain::T
end

