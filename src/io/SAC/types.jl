# Definition, constructors and other operations only to do with the SACTrace type

# Composite type for SAC evenly-spaced time series data
@eval mutable struct SACTrace
    $([:($(s)::SACFloat) for s in SAC_FLOAT_HDR]...)
    $([:($(s)::SACInt) for s in SAC_INT_HDR]...)
    $([:($(s)::SACBool) for s in SAC_BOOL_HDR]...)
    $([:($(s)::SACChar) for s in SAC_CHAR_HDR]...)
    # The time series, accessed with .t
    t::Array{SACFloat,1}
    function SACTrace(delta_in::Number, npts_in::Integer, b_in=0.)
        delta_in > 0 || error("SACTrace: delta must be positive")
        npts_in >= 0 || error("SACTrace: npts must be 0 or larger")
        # Variables are by default undefined, or false for bools
        $([:($(s) = SAC_RNULL) for s in SAC_FLOAT_HDR]...)
        $([:($(s) = SAC_INULL) for s in SAC_INT_HDR]...)
        $([:($(s) = false) for s in SAC_BOOL_HDR]...)
        $([:($(s) = SAC_CNULL) for s in SAC_CHAR_HDR]...)
        # Variables which must be present
        npts = convert(SACInt, npts_in)
        delta = convert(SACFloat, delta_in)
        b = b_in
        e = b + (npts - 1)*delta
        t = zeros(npts)
        depmin = 0.
        depmax = 0.
        nvhdr = SAC_VER_NUM
        iftype = 1
        idep = 5
        iztype = 9
        ievtyp = 5
        leven = true
        lovrok = true
        lcalda = true
        new($([:($(s)) for s in [SAC_FLOAT_HDR; SAC_INT_HDR; SAC_BOOL_HDR; SAC_CHAR_HDR]]...),
            t)
    end
end
@doc """
    SACTrace(delta, npts, b=0.) -> ::SACTrace

Construct a composite type holding an evenly-spaced SAC time-series trace, where the trace
is accessed through the field name `t`.  Supply the constant sampling interval `delta`
in seconds, and the number of points in the trace `t`.  Optionally, specify the trace
start time `b` in seconds.

    SACTrace(v::AbstractVector, delta, b=0.) -> ::SACTrace

Construct a `SACTrace` by supplying an array `v`, sampling interval `delta` and optionally
the starting time.

    SACTrace(d::Vector{UInt8}, file=""; swap=true, terse=false, check_npts=true) -> ::SACTrace

Construct a SACTrace from a raw array of bytes representing some data in SAC format.
If `swap` is false, then non-native-endian files are not converted.  If `terse` is
true, then warnings about swapping are not written.  If `check_npts` is false, then
parts of files are read without error.
""" SACTrace


function SACTrace(data::Vector{UInt8}, file=""; swap::Bool=true, terse::Bool=false,
        check_npts::Bool=true)
    len = SAC_BYTE_LEN
    clen = 2*SAC_BYTE_LEN
    # Determine endianness and act accordingly
    nvhdr = reinterpret(SACInt, data[SAC_NVHDR_POS*len+1:(SAC_NVHDR_POS+1)*len])[1]
    native = if nvhdr == SAC_VER_NUM
        true
    elseif bswap(nvhdr) == SAC_VER_NUM
        false
    else
        error("Array does not appear to be SAC data")
    end
    native && MACHINE_IS_LITTLE_ENDIAN && !swap &&
        error("Data are little-endian but `swap` is `false`.  Not attempting to swap bytes" *
            (file!="" ? " for file '$file'." : "."))
    native && MACHINE_IS_LITTLE_ENDIAN && !terse &&
        @info("Data are little-endian; byteswapping")
    byteswap(x) = native ? x : bswap(x)
    # Create an empty object
    npts_in_file = (length(data) - SAC_HEADER_LEN)÷len
    trace = SACTrace(1, npts_in_file)

    ## Read header
    # Float part
    off = 0
    for (i, field) in enumerate(SAC_FLOAT_HDR)
        setfield!(trace, field,
            byteswap(reinterpret(SACFloat, data[off+1:off+len])[1]))
        off += len
    end
    # Int part
    for (i, field) in enumerate(SAC_INT_HDR)
        setfield!(trace, field,
            byteswap(reinterpret(SACInt, data[off+1:off+len])[1]))
        off += len
    end
    # Boolean part
    for (i, field) in enumerate(SAC_BOOL_HDR)
        setfield!(trace, field,
            1 == byteswap(reinterpret(SACInt, data[off+1:off+len])[1]))
        off += len
    end
    # Character part
    # kevnm header is double length, so treat separately
    trace.kstnm = String(strip(String(data[off+1:off+clen])))
    off += clen
    trace.kevnm = String(strip(String(data[off+1:off+2clen])))
    trace.kevnm == strip(rpad(SAC_CNULL, clen)^2) && (trace.kevnm = SAC_CNULL)
    off += 2clen
    for (i, field) in enumerate(SAC_CHAR_HDR)
        i <= 2 && continue
        setfield!(trace, field,
            String(strip(String(data[1+off:clen+off]))))
        off += clen
    end

    # Check file type
    (trace.iftype in (SAC_ITIME, SAC_IXY) && trace.leven) ||
        error("only evenly-spaced time-series SAC files are supported")

    # Check length
    @assert off == SAC_HEADER_LEN
    if check_npts
        npts_in_file != trace.npts &&
            error("Number of points is not as expected: have $npts_in_file " *
                "versus npts = $(trace.npts) in header" *
                (file!="" ? " for file '$file'." : "."))
    end

    # Now read in the trace
    trace.t .= reinterpret(SACFloat, data[(SAC_HEADER_LEN+1):end])
    native || (trace.t .= bswap.(trace.t))
    update_headers!(trace)
    any(isundefined.([trace.gcarc, trace.az, trace.baz])) && update_great_circle!(trace)
    trace
end

function SACTrace(v::AbstractVector, delta::Real, b=0.0)
    s = SACTrace(delta, length(v), b)
    s.t .= v
    update_headers!(s)
    s
end

"""
    (==)(a::SACTrace, b::SACTrace) -> ::Bool

Return `true` if the traces `a` and `b` are equal (that is, have all fields the same),
and `false` otherwise.
"""
function Base.:(==)(a::SACTrace, b::SACTrace)
    for f in fieldnames(SACTrace)
        if getfield(a, f) != getfield(b, f)
            return false
        end
    end
    true
end

Base.broadcastable(x::SACTrace) = Ref(x)

"""
    getindex(A::AbstractArray{SACTrace}, s::Symbol) -> Array{typeof(A[:].s)}
    A[:s] -> Array{typeof(A[:].s)}

Return an array of values containing the header with name `s` for the SACTrace
traces.  This allows one to get all the headers values by doing A[:kstnm],
for example.
"""
Base.getindex(A::AbstractArray{SACTrace}, s::Symbol) = Array{typeof(getfield(A[1], s))}([getfield(a, s) for a in A])
Base.getindex(t::SACTrace, s::Symbol) = getfield(t, s) # Also define for single trace for consistency

"""
    setindex!(s::SACTrace, value, s::Symbol)
    s[s] = value

Set the header with name `s` (a `Symbol`) to `value`.  E.g.:

    A[:kevnm] = "Blast 1"

    setindex!(A::AbstractArray{SACTrace}, value, s::Symbol)
    A[s] = value

Set the header with name `s` for all the SACTrace traces in the array `A`.  This
allows one to set all the headers at once for a set of traces by doing e.g.:

    A[:user0] = 1:length(A)

#### Setting trace values

When setting the values of the trace, use `s[:t] = k`.  This will automatically
update headers to reflect the trace and use broadcasting internally so that
the trace cannot be set to the wrong length.  Note however that a
broadcast assignment (`s[:t] .= k`) will not update the headers; this can be
done manually with `SAC.update_headers!(s)`.

    s[:t] = 1:s[:npts] # Sets to trace to be a line with constant slope
    s[:t] = sin.(2π*time(s)) # Sets the trace to a sine function
"""
function Base.setindex!(t::SACTrace, v, s::Symbol)
    if s === :t
        t.t .= v
        update_headers!(t)
    else
        setfield!(t, s, convert(typeof(getfield(t, s)), v))
        s in (:evlo, :evla, :stlo, :stla) && update_great_circle!(t)
    end
    t
end
Base.setindex!(A::AbstractArray{SACTrace}, V, s::Symbol) = setindex!.(A, V, s)

"""
    isundefined(x) -> ::Bool

If the SAC value `x` is undefined, return `true`.

    isundefined(s::SACTrace, x::Symbol) -> :: Bool

Test whether header `x` is undefined for trace `s`.
"""
isundefined(x::SACFloat) = x == SAC_RNULL
isundefined(x::SACInt) = x == SAC_INULL
isundefined(x::SACChar) = strip(x) == strip(SAC_CNULL)
isundefined(x::SACBool) = false
isundefined(s::SACTrace, x::Symbol) = isundefined(getfield(s, x))

Base.copy(s::SACTrace) = deepcopy(s)
