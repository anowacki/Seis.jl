"Type of strings"
const _STRINGS = Dict{String,Union{String,Dict{String,String}}}

struct FileDescriptorBlock
    bigendian::Bool
    revision_number::UInt16
    trace_pointer_block_size::UInt16
    number_of_traces::UInt16
    string_terminator::String
    line_terminator::String
    trace_pointers::Vector{UInt32}
    strings::_STRINGS
end

"Determine whether we should swap bytes for this file"
function _swap(magic_number)
    if magic_number === 0x553a
        ntoh, true
    elseif magic_number === 0x3a55
        ltoh, false
    else
        error("unexpected magic two-byte number at start of file")
    end
end

function Base.read(io::IO, ::Type{FileDescriptorBlock})
    magic_number = read(io, UInt16)
    swap, bigendian = _swap(magic_number)
    revision_number = swap(read(io, UInt16))
    trace_pointer_block_size = swap(read(io, UInt16))
    number_of_traces = swap(read(io, UInt16))

    num_string_terminator_chars = read(io, UInt8)
    string_terminator = if num_string_terminator_chars == 1
        st = string(read(io, Char))
        # Unused byte
        read(io, UInt8)
        st
    elseif num_string_terminator_chars == 2
        string(read(io, Char), read(io, Char))
    else
        error("number of string terminator characters is not 1 or 2")
    end

    num_line_terminator_chars = read(io, UInt8)
    line_terminator = if num_line_terminator_chars == 1
        lt = string(read(io, Char))
        # Unused byte
        read(io, UInt8)
        lt
    elseif num_line_terminator_chars == 2
        string(read(io, Char), read(io, Char))
    else
        error("number of line terminator characters is not 1 or 2")
    end

    # Seek to trace pointer block
    seek(io, 32)

    trace_pointers = Vector{UInt32}(undef, number_of_traces)
    for i in eachindex(trace_pointers)
        trace_pointers[i] = swap(read(io, UInt32))
    end

    # Seek to the free-format string block
    seek(io, 32 + Int(trace_pointer_block_size))

    strings = _read_strings(io, swap, string_terminator, line_terminator)

    FileDescriptorBlock(
        bigendian, revision_number, trace_pointer_block_size, number_of_traces,
        string_terminator, line_terminator, trace_pointers, strings
    )
end

"Alias for the storage needed for a SEG 20-bit floating point number"
const _SEGFloat20_bitstype = NTuple{5,UInt16}

"""
    _DATA_FORMAT_CODES::Dict{UInt16,Tuple{Type,String}}

Mapping between trace descriptor data format codes and the Julia type of the
data, plus a description of that format.
"""
const _DATA_FORMAT_CODES = Dict(
    0x01 => (Int16, "16-bit fixed point"),
    0x02 => (Int32, "32-bit fixed point"),
    0x03 => (_SEGFloat20_bitstype, "20-bit floating point (SEG-D)"),
    0x04 => (Float32, "32-bit floating point (IEEE standard)"),
    0x05 => (Float64, "64-bit floating point (IEEE standard)"),
)

"""
Set of data format codes which we support in Seis (currently only IEEE floating
point types and 16- or 32-bit integer).
"""
const _SUPPORTED_DATA_FORMAT_CODES = (0x01, 0x02, 0x03, 0x04, 0x05)

"""
    TraceDescriptorBlock

Block describing a single trace.
"""
struct TraceDescriptorBlock
    bigendian::Bool
    this_block_size::UInt16
    data_block_size::UInt32
    number_of_samples::UInt32
    data_format_code::UInt8
    strings::_STRINGS
end

function Base.read(io::IO, ::Type{TraceDescriptorBlock}, fd::FileDescriptorBlock; warn=true)
    # Seek forward to the start of a 4-byte block if we're not on one
    current_pos = position(io)
    block_start_pos = if current_pos%4 == 0
        current_pos
    else
        current_pos + (4 - current_pos%4)
    end
    seek(io, block_start_pos)

    # Determine endianness of numbers
    magic_number = read(io, UInt16)
    swap, bigendian = if magic_number == 0x4422
        ltoh, false
    elseif magic_number == 0x2244
        ntoh, true
    else
        error("unexpected trace descriptor block ID (should be 0x4422 or 0x2244")
    end

    this_block_size = swap(read(io, UInt16))
    data_block_size = swap(read(io, UInt32))
    number_of_samples = swap(read(io, UInt32))
    data_format_code = read(io, UInt8)

    if !haskey(_DATA_FORMAT_CODES, data_format_code)
        error(lazy"unknown data format code: $data_format_code")
    end

    if !(data_format_code in _SUPPORTED_DATA_FORMAT_CODES)
        _, description = _DATA_FORMAT_CODES[data_format_code]
        warn && @warn "Only IEEE floating point trace data are supported; trace data are $description"
    end
    seek(io, block_start_pos + 32)

    strings = _read_strings(io, swap, fd.string_terminator, fd.line_terminator)

    # Move to the end ready to read the data
    seek(io, block_start_pos + this_block_size)

    TraceDescriptorBlock(
        bigendian, this_block_size, data_block_size, number_of_samples,
        data_format_code, strings
    )
end

"""
    SEG2Trace{T}

A single trace containing a trace descriptor block and a data block,
where the data has element type `T`.  `T` is set by the user rather than
the original data and data is converted to this element type upon reading.
"""
struct SEG2Trace{T}
    td::TraceDescriptorBlock
    data::Vector{T}
end

function Base.read(io::IO, ::Type{SEG2Trace{T}}, fd::FileDescriptorBlock; warn=true) where T
    td = read(io, TraceDescriptorBlock, fd; warn=warn)

    td.data_format_code in _SUPPORTED_DATA_FORMAT_CODES ||
        error("only IEEE floating point and integer data types are supported")
    traceT, _ = _DATA_FORMAT_CODES[td.data_format_code]
    # This is a float
    bytes_per_sample = if traceT == _SEGFloat20_bitstype
        # Must have a multiple of four samples in this format
        if td.number_of_samples%4 != 0
            error("20-bit sample format requires the number of samples be divisible by four")
        end
        # 2.5 bytes per sample because four two-byte samples take up 10 bytes
        2.5
    else
        sizeof(traceT)*1.0
    end

    if (sizeof(T) < bytes_per_sample) && traceT != NTuple{5,UInt16}
        warn && @warn "chosen trace eltype ($T) is smaller than file trace eltype ($traceT); precision will be lost"
    end

    # Should exactly convert to an integer
    bytes_to_read = Int(td.number_of_samples*bytes_per_sample)
    if bytes_to_read != td.data_block_size
        error("data block size does not match number of samples and byte length of samples")
    end

    bytes = read(io, bytes_to_read)
    sample_data = reinterpret(traceT, bytes)
    if td.bigendian
        sample_data .= ntoh.(sample_data)
    else
        sample_data .= ltoh.(sample_data)
    end

    data = if traceT == _SEGFloat20_bitstype
        _data = Vector{T}(undef, td.number_of_samples) .|> bswap
        for i in 1:(td.number_of_samples÷4)
            segfloat20s = SEGFloat20(sample_data[i])
            isample1 = 4*(i - 1) + 1
            _data[isample1], _data[isample1 + 1], _data[isample1 + 2], _data[isample1 + 3] =
                _segfloat20_samples(sample_data[i], T)
        end
        _data
    else
        sample_data
    end

    SEG2Trace{T}(td, data)
end

"""
    SEG2File{T}

A single SEG2 file containing the overall file descriptor block,
and multiple [`SEG2Trace`](@ref)s.  Each one must have the same element
type `T`.
"""
struct SEG2File{T}
    fd::FileDescriptorBlock
    traces::Vector{SEG2Trace{T}}
end

function Base.read(io::IO, ::Type{SEG2File{T}}; warn=true) where T
    fd = read(io, FileDescriptorBlock)

    traces = Vector{SEG2Trace{T}}(undef, fd.number_of_traces)
    for i in eachindex(traces, fd.trace_pointers)
        seek(io, fd.trace_pointers[i])
        traces[i] = read(io, SEG2Trace{T}, fd; warn=warn)
    end

    SEG2File{T}(fd, traces)
end

function Base.show(io::IO, mime::MIME"text/plain", f::SEG2File{T}) where T
    println(io, "SEG2File{$T} ($(f.fd.number_of_traces) traces):")
    show(io, mime, f.fd)
end


"""
    _read_strings(io, swap, string_terminator, line_terminator) -> $(_STRINGS)

Read the free-format stings from the end of file and trace descriptor
blocks.  `swap` is a function which correctly byte-swaps the string
offset `UInt16` values, `string_terminator` is the single- or
double-character string which terminates each string, and `line_terminator`
is the single- or double-character string ending a line within a string.
"""
function _read_strings(io, swap, string_terminator, line_terminator)
    strings = _STRINGS()
    next_string_offset = swap(read(io, UInt16))

    while next_string_offset > 0
        _string = readuntil(io, string_terminator)
        key_end_index = findfirst(in((' ', Char(20))), _string)
        key = _string[1:key_end_index-1]

        if key == "NOTE"
            # Special case of notes which have several key-value pairs
            notes = Dict{String,String}()

            for note_string in split(_string[key_end_index+1:end], line_terminator)
                isempty(strip(note_string)) && continue
                note_key_value = match(r"^\s*([A-Z_]+)\s+(.*)\s*$", note_string)
                if isnothing(note_key_value) || length(note_key_value.captures) != 2
                    strings[key] = String(note_string)
                else
                    note_key, note_value = note_key_value
                    notes[note_key] = strip(note_value)
                end
            end

            strings[key] = notes
        else
            value = _string[key_end_index+1:end]
            strings[key] = value
        end

        next_string_offset = swap(read(io, UInt16))
    end

    strings
end

function Base.show(io::IO, ::MIME"text/plain", o::T) where {T<:Union{FileDescriptorBlock,TraceDescriptorBlock}}
    print(io, T)
    for f in fieldnames(T)
        val = getfield(o, f)
        print(io, '\n')
        if val isa String
            print(io, " ", f, ": ", repr(val))
        else
            print(io, " ", f, ": ", (val isa Unsigned ? Int(val) : val))
        end
    end
end
