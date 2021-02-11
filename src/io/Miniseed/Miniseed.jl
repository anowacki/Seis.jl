"""
# `Miniseed`

Module for dealing with miniseed data.
"""
module Miniseed

import SeisIO
import ..Seis
using ..Seis: SeisIOIO, Trace, CartTrace, AbstractTrace

# Default miniseed IO parameters used by SeisIO
const NX_NEW = SeisIO.KW.nx_new
const NX_ADD = SeisIO.KW.nx_add

DEFAULT_TRACE = Trace{Float32, Vector{Float32}, Seis.Geographic{Float32}}

"""
    read(file; maximum_gap, maximum_total_offset, kwargs...) -> ::Vector{$DEFAULT_TRACE}
    read(T, file; kwargs...) -> ::Vector{T}

Read miniSEED data from `file` as a set of `Trace`s, optionally specifying the
type of `Trace` `T <: AbstractTrace`.

For details of keywords arguments, see [`Seis.read_mseed`](@ref).
"""
function read(::Type{T}, file;
        maximum_gap=nothing, maximum_offset=nothing,
        nx_new=NX_NEW, nx_add=NX_ADD, memmap=false, strict=false, verbose=0
        ) where {T<:AbstractTrace}
    seisdata = SeisIO.SeisData()
    SeisIO.SEED.read_mseed_file!(seisdata, file, Int64(nx_new), Int64(nx_add),
        memmap, strict, verbose)
    SeisIOIO.parse_seisio(T, seisdata, file; maximum_gap=maximum_gap,
        maximum_offset=maximum_offset)
end

"""
    read(data::AbstractVector{UInt8}; kwargs...) -> ::Vector{$DEFAULT_TRACE}
    read(T, data; kwargs...) -> ::Vector{T}

Read miniSEED data from memory as a set of `Trace`s, optionally specifying the
type of `Trace` `T <: AbstractTrace`.

For details of keywords arguments, see [`Seis.read_mseed`](@ref).
"""
function read(::Type{T}, data::AbstractVector{UInt8};
        maximum_gap=nothing, maximum_offset=nothing,
        nx_new=NX_NEW, nx_add=NX_ADD, strict=false, verbose=0
        ) where {T<:AbstractTrace}
    seisdata = SeisIO.SeisData()
    io = IOBuffer(data)
    SeisIO.SEED.parsemseed!(seisdata, io, Int64(nx_add), Int64(nx_new),
        strict, verbose)
    SeisIOIO.parse_seisio(T, seisdata; maximum_gap=maximum_gap,
        maximum_offset=maximum_offset)
end

read(file; kwargs...) = read(DEFAULT_TRACE, file; kwargs...)

end # module
