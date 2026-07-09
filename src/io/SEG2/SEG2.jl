"""
# `SEG2`

Module for reading and writing files in the SEG2 format.
"""
module SEG2

import Dates
import ..Seis

include("float20.jl")
include("types.jl")
include("util.jl")

"""
    read_seg2(io_or_file[, T=Float32]; warn=true) -> seg2file::SEG2.SEG2File

Read the data and headers from a SEG2 file, setting the trace data element
type to `T`.  Note that this will convert trace data in a different format.
"""
function read_seg2(io::IO, ::Type{T}; warn=true) where T
    read(io, SEG2File{T}; warn=warn)
end
read_seg2(io::IO; kwargs...) = read_seg2(io, Float32; kwargs...)

read_seg2(file, ::Type{T}; kwargs...) where T = open(io -> read_seg2(io, T; kwargs...), file)
read_seg2(file; kwargs...) = read_seg2(file, Float32; kwargs...)

end # module
