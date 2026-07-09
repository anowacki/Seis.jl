"""
    SEGFloat20

Struct holding 10 bytes which represent four values in SEG 20-bit
floating point format.

These can be converted to other floating point formats 
"""
struct SEGFloat20
    bytes::NTuple{5,UInt16}
end

"""
    _segfloat20_samples(x::SEGFloat20, ::Type{T}=Int32) where T -> (sample1, sample2, sample3, sample4)
    _segfloat20_samples(x::NTuple{5,UInt16}, ::Type{T}=Int32) where T -> (sample1, sample2, sample3, sample4)

Convert a set of five `UInt16`s storing four samples in SEG 20-bit floating
point format into separate samples of type `T`.
"""
function _segfloat20_samples(x::NTuple{5,UInt16}, ::Type{T}=Int32) where T
    # Exponents for each of the four samples
    exponents = _exponents(x[1])

    ntuple(4) do i
        T(_segfloat20_value(x[i + 1], exponents[i]))
    end
end
function _segfloat20_samples(x::SEGFloat20, ::Type{T}=Float32) where T
    _segfloat20_samples(x.bytes, T)
end

"Return four exponents from a single 16-bit unsigned int"
function _exponents(x::UInt16)
    ntuple(4) do i
        (x & (0x000f << (4*(i - 1)))) >> (4*(i - 1))
    end
end

"Combine a 4-bit exponent and 16-bit signed int to give the floating-point value"
function _segfloat20_value(byte::UInt16, exponent::UInt16)
    sign = (byte & 0b1000_0000_0000_0000) > 0
    value = Int32(reinterpret(Int16, byte + sign))*(Int32(1) << exponent)
end
