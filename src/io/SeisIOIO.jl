"""
    SeisIOIO

Module for common tools for interacting with SeisIO types and functions.

The name of the module conveys that it deals with IO via SeisIO.
"""
module SeisIOIO

using Dates: unix2datetime

import SeisIO
import ..Seis

"""
    Gap

Struct representing a gap or overlap of data in a SeisIO.SeisChannel.
"""
struct Gap
    "Sample of start of next piece of continuous data"
    start_sample::Int
    "Offset from previous sample in μs"
    offset_μs::Float64
end

"Find gaps in SeisIO.SeisChannels"
gaps(chan::SeisIO.SeisChannel) = (Gap(chan.t[i,1], chan.t[i,2]) for i in 2:(size(chan.t, 1)-1))

"Sum total of offsets in s."
total_offset(chan::SeisIO.SeisChannel) = sum(x->x.offset_μs, gaps(chan))/1_000_000

"Largest absolute offset in s"
largest_offset(chan) = maximum(x->abs(x.offset_μs), gaps(chan))/1_000_000

"""
    parse_seisio(T, seisdata, file=missing; maximum_gap, maximum_offset) -> ::Vector{T}
    parse_seisio(seisdata, file=missing; maximum_gap, maximum_offset) -> ::Vector{Trace{Float64}}

Transform data in `SeisIO.SeisData` form into `Trace`s for use with Seis.

Specify the type of trace `T`, which defaults to a `Float64` geographic Trace.

See note in [`Seis.read_mseed`](@ref) on handling of gapped/overalapped data.

If `maximum_gap` is specified, then traces are only split at gaps when
the absolute gap length in s is greater than `maximum_gap`.
If `maximum_offset` is specified, traces are split when the sum total of
offsets (i.e., total of negative and positive gaps) in s is greater than
`maximum_offset`.
"""
function parse_seisio(::Type{T}, seisdata::SeisIO.SeisData, file=missing;
        maximum_gap=nothing, maximum_offset=nothing) where {T<:Seis.AbstractTrace}
    nchannels = length(seisdata)
    traces = T[]
    for i in 1:nchannels
        channel = seisdata[i]
        delta = 1/channel.fs
        time = seisio_date(channel)
        # Deal with gaps
        max_gap = maximum_gap === nothing ? delta : maximum_gap
        max_offset = maximum_offset === nothing ? delta : maximum_offset
        if SeisIO.ngaps(channel.t) == 0 ||
                (total_offset(channel) <= max_offset && largest_offset(channel) <= max_gap)
            t = T(0, delta, channel.x)
            assign_name!(t, channel)
            assign_position!(t, channel)
            t.evt.time = time
            push!(traces, t)
        else
            data, bs = chunks(channel)
            ts = [T(b, delta, d) for (d, b) in zip(data, bs)]
            assign_name!.(ts, (channel,))
            assign_position!.(ts, (channel,))
            ts.evt.time = time
            append!(traces, ts)
        end
    end
    traces.meta.mseed_file = file
    traces
end

parse_seisio(seisdata::SeisIO.SeisData, file=missing; kwargs...) =
    parse_seisio(Seis.Trace{Float32,Vector{Float32},Seis.Geographic{Float32}}, seisdata, file; kwargs...)
parse_seisio(T::Type, chan::SeisIO.SeisChannel, file=missing; kwargs...) =
    parse_seisio(T, SeisIO.SeisData(chan), file; kwargs...)
parse_seisio(chan::SeisIO.SeisChannel, file=missing; kwargs...) =
    parse_seisio(Seis.Trace{Float32, Vector{Float32}, Seis.Geographic{Float32}}, chan, file; kwargs...)

"Convert the id field of a SeisIO.SeisChannel into the name fields of a Trace."
function assign_name!(t, channel::SeisIO.SeisChannel)
    tokens = split(channel.id, '.')
    if length(tokens) == 4
        t.sta.net, t.sta.sta, t.sta.loc, t.sta.cha = tokens
    else
        @warn("channel code is not in an expected format")
        t.sta.sta = channel.id
    end
    t
end

assign_position!(t::Seis.AbstractTrace, x) = assign_position!(t.sta, x)

"Convert the loc field of a SeisIO.SeisChannel into the position field
of a Trace's Station"
assign_position!(sta::Seis.Station, channel::SeisIO.SeisChannel) = assign_position!(sta, channel.loc)

"Fall back method doesn't assign anything"
assign_position!(sta::Seis.Station, loc) = nothing

"Assign geographic locations to geographic stations"
function assign_position!(sta::Seis.GeogStation, loc::SeisIO.GeoLoc)
    # Can only assume that everything being 0 means these values aren't set
    if all(x -> x==0, (loc.lat, loc.lon, loc.el, loc.dep, loc.az, loc.inc))
        return nothing
    end
    # Otherwise, at least one field is not 0 and we assume the others are set
    sta.lon = loc.lon
    sta.lat = loc.lat
    sta.elev = loc.el
    sta.dep = loc.dep/1000
    sta.azi = loc.az
    sta.inc = loc.inc
    loc.datum != "WGS84" && (sta.meta.datum = loc.datum)
    nothing
end

"Assign Cartesian locations to Cartesian stations"
function assign_position!(sta::Seis.CartStation, loc::SeisIO.XYLoc)
    if all(x -> x==0, (loc.x, loc.y, loc.z, loc.az, loc.inc))
        return nothing
    end
    sta.x = loc.x
    sta.y = loc.y
    sta.z = loc.z
    sta.azi = loc.az
    sta.inc = loc.inc
    sta.meta.origin = (x=loc.ox, y=loc.oy, z=loc.oz)
    sta.meta.datum = loc.datum
    nothing
end

"""
    seisio_date(chan:SeisIO.SeisChannel) -> ::Dates.DateTime

Find the date of the start of a SeisIO.SeisChannel.
"""
seisio_date(chan::SeisIO.SeisChannel) = unix2datetime(chan.t[1,2]/1e6)

"""
    chunks(chan::SeisIO.SeisChannel) -> data, begin_time

Divide the data in `chan` into vectors of `data`, each with an associated
`begin_time` which is relative to the channel's start time.
"""
function chunks(channel::SeisIO.SeisChannel)
    delta = 1/channel.fs
    ngaps = SeisIO.ngaps(channel.t)
    if ngaps == 0
        return [channel.x], [0.0]
    else
        data = Vector{Vector{Float64}}(undef, ngaps + 1)
        bs = Vector{Float64}(undef, ngaps + 1)
        for i in 1:(ngaps+1)
            start_sample = channel.t[i,1]
            end_sample = channel.t[i+1,1] - 1
            if i == ngaps + 1
                end_sample += 1
            end
            offset_μs = i == 1 ? 0 : channel.t[i,2]
            data[i] = channel.x[start_sample:end_sample]
            bs[i] = (start_sample - 1)*delta + offset_μs/1e6
        end
        return data, bs
    end
end

end # module
