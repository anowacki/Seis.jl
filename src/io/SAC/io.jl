# Input/output routines

"""
    read(file; swap=true, terse=false, check_npts=true, header_only=false) -> s::SACTrace

Return the SAC trace as read from file `file` as a `SACTrace` object.  If `swap` is false,
then auto-byteswapping is not performed and an error is returned if the file is not
of the assumed endianness.  Autoswapping is reported unless `terse` is `true`.

By default an error is thrown if the trace on disk is not the same length
as specified by the `npts` header variable.  If `check_npts` is `false`,
then no error is thrown and the trace is returned.  (Note that in this
case the `npts` header is out of sync with the trace.)  This option may be
useful for SAC data which has been incorrectly written.

If `header_only` is true, then only the header part of the file is read
from disk.  This option sets `check_npts` to `false`.
"""
function read(file; swap::Bool=true, terse::Bool=false, header_only::Bool=false,
        check_npts::Bool=!header_only)
    open(file, "r") do io
        read(io; swap=swap, terse=terse, header_only=header_only,
            check_npts=check_npts, file=file)
    end
end

function read(io::IO; swap::Bool=true, terse::Bool=true, header_only::Bool=false,
        check_npts::Bool=!header_only, file="")
    data = header_only ? Base.read(io, SAC_HEADER_LEN) : Base.read(io)
    header_only && (check_npts = false)
    SACTrace(data, file; swap=swap, terse=terse, check_npts=check_npts)
end

"""
    read_cut(file, b, e; swap=true, terse=false) -> s::SACTrace
    read_cut(file, b_header, b_time, e_header, e_time; swap=true, terse=false) -> s::SACTrace

Return the trace `s`, cut between beginning time `b` and end time `e` in seconds
relative to the zero time.

Optionally specify headers to which `b_time` and `e_time` are relative as
`b_header` and `e_header` respectively.  For example:

```
s = read_cut(file, :a, -1, :f + 20)
```

If `b` is before the start of the trace, then the trace is read from the beginning;
similarly, if `e` is after the end of the trace, the trace is read up until the end.
In either case, a warning is issued unless `terse` is `true`.
"""
function read_cut(io::IO, hb::Symbol, b::Real, he::Symbol, e::Real; swap=true, terse=false, file=nothing)
    b < e || throw(ArgumentError("Start cut must be before stop cut"))
    len = SAC_BYTE_LEN
    clen = 2*SAC_BYTE_LEN
    header = Base.read(io, SAC_HEADER_LEN)
    # Determine endianness and act accordingly
    nvhdr = reinterpret(SACInt, header[SAC_NVHDR_POS*len+1:(SAC_NVHDR_POS+1)*len])[1]
    native = if nvhdr == SAC_VER_NUM
        true
    elseif bswap(nvhdr) == SAC_VER_NUM
        false
    else
        error("File does not appear to be SAC data (nvhdr = $nvhdr)")
    end
    native || swap || error("File $file is non-native-endian but `swap` is `false`")
    # Create new SACTrace object, setting the trace to have one zero value
    append!(header, zeros(UInt8, 4))
    s = SACTrace(header, check_npts=false)
    # Check args
    isundefined(s, hb) && error("Header $hb is undefined")
    isundefined(s, he) && error("Header $he is undefined")
    b < s.b && !terse &&
        @warn("Start cut $b is before trace start $(s.b); reading from trace start")
    e > s.e && !terse &&
        @warn("End cut $e is after trace end $(s.e); reading to trace end")
    # Now do cutting
    ib = round(Int, (b - getfield(s, hb))/s.delta) + 1
    ie = round(Int, (e - getfield(s, he))/s.delta) + 1
    npts = ie - ib + 1
    @assert ib >= 1 && ie <= s.npts
    s.b += (ib - 1)*s.delta
    s.npts = npts
    seek(io, SAC_HEADER_LEN + len*(ib - 1))
    s.t = reinterpret(SACFloat, Base.read(io, len*npts))
    native || (s.t .= bswap.(s.t))
    update_headers!(s)
end
read_cut(file, b::Real, e::Real; kwargs...) = read_cut(file, :b, b, :b, e; kwargs...)

read_cut(file, hb, b, he, e; swap=true, terse=false) = open(file, "r") do io
    read_cut(io, hb, b, he, e; swap=swap, terse=terse, file=file)
end

"""
    file_is_native_endian(file)

Return `true` if `file` is a native-endian (defined by a constant in the module to be
little-endian on this machine) SAC file, and `false` if not.

The heuristic is thus: native-endian files have bytes 305:308 which are
a representation of the value `6`.  `6` is the current magic SAC file version number,
hard-coded into the routine.
"""
function file_is_native_endian(file::String)
    nvhdr = try
        d = open(file, "r") do f
            seek(f, SAC_NVHDR_POS*SAC_BYTE_LEN)
            reinterpret(SACInt, Base.read(f, SAC_BYTE_LEN))[1]
        end
    catch err
        error("SAC.file_is_native_endian: Cannot open file '$file' for reading " *
            "or file is not the correct type (error $err)")
    end
    if nvhdr == SAC_VER_NUM
        return true
    elseif bswap(nvhdr) == SAC_VER_NUM
        return false
    else
        error("SAC.file_is_native_endian: File '$file' does not appear to be a " *
            "valid SAC file (nvhdr is $nvhdr)")
    end
end

"Write floats or integers either swapped or native-endian"
_write_swap(swap::Bool, F::IO, x) = Base.write(F, swap ? Base.bswap.(x) : x)
"Write a String as a SAC string of the correct length, padded with ' 's"
_write_string(F::IO, x::String, maxlen::Integer) =
    Base.write(F, x[1:min(length(x),maxlen)]*" "^(max(0,maxlen-length(x))))

# Write the header part of a `SACTrace` to `io`.  If `byteswap` is `true`
# then the header is written littleendian.
@eval function _write_header(s::SACTrace, io; byteswap=SAC_FORCE_SWAP)
    # Write header
    $([:(_write_swap(byteswap, io, s.$s)) for s in [SAC_FLOAT_HDR; SAC_INT_HDR]]...)
    $([:(_write_swap(byteswap, io, SACInt(s.$s))) for s in SAC_BOOL_HDR]...)
    # No byte-swapping needed for characters, but pad them to the correct length
    _write_string(io, s.kstnm, SACCHARLEN)
    # Handle special case of double-length kevnm  header
    kevnm = isundefined(s.kevnm) ? "-12345  -12345  " : s.kevnm
    _write_string(io, kevnm, 2*SACCHARLEN)
    $([:(_write_string(io, s.$s, SACCHARLEN)) for s in SAC_CHAR_HDR[3:end]]...)
end

"""
    write(s::SACTrace, file; byteswap)
    write(S::AbstractArray{SACTrace}, files; byteswap)

Write a SAC trace `s` to `file`, or a set of traces `S` to a set of files `files`.
Set `byteswap` to `false` to force writing in native-endian format; set to `true`
to write bigendian files (MacSAC type).  The default is to write bigendian format.
"""
function write(s::SACTrace, file; byteswap=SAC_FORCE_SWAP)
    open(file, "w") do io
        write(s, io; byteswap=byteswap)
    end
end

"""
    write(s::SACTrace, io::IO; byteswap)

Write a SAC file to `io`, which may be an `IOStream` (such as a file handle from
opening a file), `IOBuffer` or similar.
"""
function write(s::SACTrace, io::IO; byteswap=SAC_FORCE_SWAP)
    _write_header(s, io; byteswap=byteswap)
    _write_swap(byteswap, io, s.t)
end


function write(s::AbstractArray{SACTrace}, file::Array{String}; args...)
    length(s) == length(file) || error("SAC.write: Arrays must be same length")
    for i = 1:length(s)
        write(s[i], file[i]; args...)
    end
    return
end

"""
    overwrite_header(s::SACTrace, file; check=true, byteswap=$(SAC_FORCE_SWAP))

For an existing SAC `file` on disk, overwrite its header with the
header of the trace `s`.

When `check` is `true` (the default), `overwrite_header` checks that
`file` is an existing SAC trace with the correct number of points
in the data trace, and will determine the file endianness in order
to write the headers correctly.  However, if `check` is `false`,
then no checks are made.  In this case, the header will be written
in bigendian endianness (MacSAC or SAC/BRIS format) unless
`byteswap` is `true`.  `byteswap` has no effect if `check` is `true`.

This function is useful for updating headers for files on disk
without having to read and write the entire trace from and to
the disk.

!!! warning
    With the option `check=false`, no check is made that the header
    written to disk matches the trace on disk in any way.  Take
    care in particular to write a value of NPTS to the header
    which matches the number of points in the pre-existing SAC file.

    `overwrite_header` will happily overwrite the first bytes of **ANY**
    file you point it at if `check` is `false`, making no check that
    `file` is actually a SAC file.

!!! note
    The NPTS header is only taken from the field in the header of `s`.
    The length of the trace in memory is ignored.
"""
function overwrite_header(s::SACTrace, file::AbstractString; check=true, byteswap=SAC_FORCE_SWAP)
    if check
        isfile(file) || error("file '$file' does not exist")
        npts = (filesize(file) - SAC_HEADER_LEN)÷4

        if npts != s.npts
            error("trace npts and file npts do not match")
        end

        byteswap = !file_is_native_endian(file)
    end
        
    open(file, "a") do io
        seekstart(io)
        _write_header(s, io; byteswap=byteswap)
    end

    return
end

"""
    read_wild(pat, dir="."; echo=true, header_only=false) -> sac_traces, files

Read files matching globbing pattern `pat` from directory `dir`.
If `echo` is false, do not show which files are being read.
Read only trace headers is `header_only` is `true`.

Returns an array of `SACTrace`s `sac_traces`, and an array of file names `files`.
**NB:** If no files are matched, then empty `SACTrace` and `String` arrays
are returned and a warning message printed.
"""
function read_wild(pat::String, dir::String=".";
        echo::Bool=true, header_only::Bool=false)
    files = Glob.glob(pat, dir)
    n = size(files, 1)
    if n == 0
        @info("SAC.read_wild: No files matching '$pat' in directory '$dir'")
        return SACTrace[], String[]
    end
    A = Array{SACTrace}(undef, n)
    for i = 1:n
        echo && @info("SAC.read: '$(files[i])'")
        A[i] = SAC.read(files[i]; terse=!echo, header_only=header_only)
    end
    return A, files
end

# List printing
function Base.show(io::IO, s::SACTrace)
    out = "SAC.SACTrace(delta=$(s[:delta])"
    for hdr in (:b, :npts, :kstnm, :gcarc, :az, :baz)
        if !isundefined(s[hdr])
            out *= ", " * string(hdr) * "=" * strip(string(s[hdr]))
        end
    end
    print(io, out * ")")
end

# Printing to the REPL, for instance
function Base.show(io::IO, ::MIME"text/plain", s::SACTrace)
    hdr_string_len = maximum(length.(string.(SAC_ALL_HDR)))
    print(io, "SAC.SACTrace:")
    for hdr in SAC_ALL_HDR
        !isundefined(s[hdr]) &&
            print(io, "\n", lpad(string(hdr), hdr_string_len, ' ') * ": ", s[hdr])
    end
end
