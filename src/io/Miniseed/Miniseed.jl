"""
# `Miniseed`

Module for dealing with miniseed data.
"""
module Miniseed

import Dates
import LibMseed
using NanoDates: NanoDates, NanoDate
import ..Seis
using ..Seis: Trace, CartTrace, AbstractTrace

const DEFAULT_TRACE = Trace{Float64, Vector{Float32}, Seis.Geographic{Float64}}

"""
    read(file; maximum_gap=nothing, verbose=0) -> ::Vector{$DEFAULT_TRACE}
    read(file, T; kwargs...) -> ::Vector{T}

Read miniSEED data from `file` as a set of `Trace`s, optionally specifying the
type of `Trace` `T <: AbstractTrace`.

For details of keywords arguments, see [`Seis.read_mseed`](@ref).
"""
function read(file, ::Type{T};
    header_only=false,
    maximum_gap=nothing,
    verbose=0,
) where {T<:AbstractTrace}
    tracelist = LibMseed.read_file(file;
        headers_only=header_only,
        time_tolerance=maximum_gap,
        verbose_level=verbose,
    )
    parse_tracelist(T, tracelist, file; header_only=header_only)
end

"""
    read(data::AbstractVector{UInt8}; kwargs...) -> ::Vector{$DEFAULT_TRACE}
    read(T, data; kwargs...) -> ::Vector{T}

Read miniSEED data from memory as a set of `Trace`s, optionally specifying the
type of `Trace` `T <: AbstractTrace`.

For details of keywords arguments, see [`Seis.read_mseed`](@ref).
"""
function read(data::AbstractVector{UInt8}, ::Type{T};
        header_only=false,
        maximum_gap=nothing,
        verbose=0,
) where {T<:AbstractTrace}
    tracelist = LibMseed.read_buffer(data;
        headers_only=header_only,
        time_tolerance=maximum_gap,
        verbose_level=verbose,
    )
    parse_tracelist(T, tracelist; header_only=header_only)
end

read(file; kwargs...) = read(file, DEFAULT_TRACE; kwargs...)

"""
    write(file, trace; verbose=0, pubversion=1, record_length=nothing, version=2)

Write the data contained in the `AbstractTrace` `trace` to `file` on disk
in miniSEED format.

If `trace` does not have an origin time set (in the `.evt.time` field),
an error is thrown.

For keyword arguments, see [`Seis.write_mseed`](@ref).

!!! note
    Seis uses `NanoDates.NanoDate`s to represent points in time, giving nanosecond
    resolution.  miniSEED version 2 uses microseconds, therefore
    trace start times will be truncated to the nearest microsecond when
    writing to this format.  It is currently the default because although
    miniSEED version 3 uses nanoseconds, is not yet widely used.
"""
function write(file, t::AbstractTrace;
        append=false, verbose=0, pubversion=1, record_length=nothing, version=2)

    ismissing(t.evt.time) &&
        throw(ArgumentError("trace has no origin time set; cannot write miniSEED"))
    b = Seis.starttime(t)
    # Start offset to the nearest nanosecond
    ns = Dates.Nanosecond(round(Int64, b*1_000_000_000, RoundToZero))
    # First sample date to nanosecond precision
    startnanodate = Seis.origin_time(t) + ns
    # Nanosecond precision of first sample
    startnanotime = _nanodate2libmseed_nanodatetime(startnanodate)
    id = trace_id(t.sta)
    sample_rate = inv(t.delta)
    LibMseed.write_file(file, Seis.trace(t), sample_rate, startnanotime, id;
        append=append, verbose_level=verbose, pubversion=pubversion,
        record_length=record_length)
end

"""
    write(file, traces; kwargs...)

Write multiple `traces` to a single miniSEED file.
"""
function write(file, traces::AbstractArray{<:AbstractTrace}; kwargs...)
    isempty(traces) && return
    for (i, t) in enumerate(traces)
        append = i > 1
        write(file, t; kwargs..., append=append)
    end
end

"""
    parse_tracelist(T, tracelist::LibMseed.MseedTraceList, file=nothing; header_only=false) -> ::Vector{T}

Convert a `LibMseed.MseedTraceList` into a `Vector{T}`, where `T <: AbstractTrace`.

The channel ID of each channel is converted to standard SEED network,
station, location and channel codes as best as possible, with all of these
being potentially empty.

!!! note
    Seis uses `NanoDates.NanoDate`s to represent points in time, giving nanosecond
    resolution.  miniSEED version 2 uses microseconds, while version 3
    uses nanoseconds, though the newer version is not yet widely used.
"""
function parse_tracelist(
    ::Type{T},
    tracelist::LibMseed.MseedTraceList,
    file=nothing;
    header_only=false
) where {T<:AbstractTrace}
    ntraces = sum(x->length(x.segments), tracelist.traces)
    traces = Vector{T}(undef, ntraces)
    itrace = 0

    for mstrace in tracelist.traces
        id = LibMseed.channel_code_parts(mstrace)
        # Traces with no channel code are sometimes written by other software
        # so that their code is all spaces.
        cha = all(isspace, id.cha) ? missing : id.cha
        # We don't expect empty codes for anything apart from loc,
        # so if any are empty, assume empty fields mean there is no code
        if any(isempty, (id.net, id.sta)) || cha === missing
            loc = isempty(id.loc) ? missing : id.loc
        # Otherwise allow loc to be empty but not missing
        else  
            loc = id.loc
        end
        net = isempty(id.net) ? missing : id.net
        sta = isempty(id.sta) ? missing : id.sta

        for segment in mstrace.segments
            itrace += 1
            starttime = _libmseed_nanodatetime2nanodate(segment.starttime)
            delta = 1/segment.sample_rate
            tr = T(0, delta, segment.data)
            tr.sta.net = net
            tr.sta.sta = sta
            tr.sta.loc = loc
            tr.sta.cha = cha
            tr.evt.time = starttime
            if file !== nothing
                tr.meta.mseed_file = file
            end

            if header_only
                tr.meta.mseed_nsamples = segment.sample_count
                tr.meta.mseed_enddate = _libmseed_nanodatetime2nanodate(segment.endtime)
                tr.meta.mseed_endtime = delta*(segment.sample_count - 1)
            end

            traces[itrace] = tr
        end

    end
    traces
end

"""
    trace_id(sta::Station, maxlengths=(2, 5, 2, 3)) -> id::String
    trace_id(net, sta, cha, loc, maxlengths=(2, 5, 2, 3)) -> id::String

Convert the network, station, location and channel code from a `Seis.Station`
or as individual strings into a single string `id` for use when saving
a channel in miniSEED format.
"""
function trace_id(sta::Seis.Station, maxlengths=(2, 5, 2, 3))
    length(maxlengths) == 4 && all(x -> x>=0, maxlengths) ||
        throw(ArgumentError("maxlengths should have four elements all 0 or more"))
    net, sta, loc, cha = check_id(sta, maxlengths)
    # `check_id` checks that strings are ASCII so we are safe with indexing below
    net = first(net, maxlengths[1])
    sta = first(sta, maxlengths[2])
    loc = first(loc, maxlengths[3])
    # If channel code is fewer than 3 characters, pad it with spaces;
    # libmseed is happy with empty bits for all the other id parts but
    # likes spaces for unused band, source and position fields.
    cha = lpad(first(cha, maxlengths[4]), maxlengths[4], ' ')
    "FDSN:" * join((net, sta, loc, cha...), '_')
end

"""
    check_id(net, sta, loc, cha, maxlengths=(2, 5, 2, 3)) -> net, sta, loc, cha
    check_id(sta, maxlengths=(2, 5, 2, 3)) -> net, sta, loc, cha

Check that the `net`, `sta`, `loc` and `cha` fields from a `Seis.Station`
are the correct lengths, warning if any are too long or contain non-ASCII
characters.
"""
check_id(sta::Seis.Station, maxlengths) =
    check_id(sta.net, sta.sta, sta.loc, sta.cha, maxlengths)

function check_id(net, sta, loc, cha, maxlengths)
    parts = (net, sta, loc, cha) = coalesce.((net, sta, loc, cha), "")
    names = ("network", "station", "location", "channel")
    for part in parts
        isascii(part) ||
            throw(ArgumentError("string \"$part\" in station ID is not " *
                                "ASCII, and cannot be written as miniSEED"))
    end
    # Check for id elements which are too long for this version of miniSEED
    for (part, maxlength, name) in zip(parts, maxlengths, names)
        if length(part) > maxlength
            @warn("$name code truncated to first $maxlength characters")
        end
    end

    parts
end

"""
    _libmseed_nanodatetime2nanodate(dt::LibMseed.NanosecondDateTime) -> nd::NanoDates.NanoDate

Convert a `LibMseed.NanosecondDateTime` into a `NanoDates.NanoDate`.
"""
function _libmseed_nanodatetime2nanodate(dt::LibMseed.NanosecondDateTime)
    NanoDate(LibMseed.datetime(dt)) + LibMseed.nanoseconds(dt)
end

"""
    _nanodate2libmseed_nanodatetime(nd::NanoDate) -> nst::LibMseed.NanosecondDateTime

Convert a `NanoDates.NanoDate` into a `LibMseed.NanosecondDateTime
"""
function _nanodate2libmseed_nanodatetime(nd::NanoDate)
    LibMseed.NanosecondDateTime(
        Dates.DateTime(nd),
        Dates.Nanosecond(1000*Dates.microsecond(nd) + Dates.nanosecond(nd))
    )
end

end # module
