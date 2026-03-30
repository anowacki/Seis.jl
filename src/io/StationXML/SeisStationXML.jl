"""
# `SeisStationXML`

Module for reading and writing StationXML files into Seis structures
"""
module SeisStationXML

import StationXML
import ..Seis

"""
    read(file, S::Type{Seis.Station}; warn=true, full=false) -> stations::Vector{S}
    read(io, S::Type{Seis.Station}; warn=true, full=false) -> stations::Vector{S}

Read a StationXML file from `file` on disk or a device `io::IO`, optionally
specifying a type `S` which must be a `Station` type.
By default, `S` is `GeogStation{Float64}`.
Return a vector of all channels specified in the StationXML file.

Pass `warn=false` to turn off warnings when encountering unexpected fields
(which are ignored) in the StationXML file.

Several fields in each station's `.meta` field are filled.  If `full` is
`true` (the default), the station's `meta.stationxml` field contains a full
`StationXML.FDSNStationXML` object for that channel.

!!! note
    For channels with full response information, the `FDSNStationXML` field
    in `meta.stationxml` may be large, which is why `full=false` by default.

# Fields stored
Upon return, each of the `stations` has its fields filled as follows from
the StationXML file information:
- `meta.burial_depth`: Station burial depth below the surface in m
- `meta.startdate`: Date at which channel started recording
- `meta.enddate`: Date at which channel stopped recording, if any
- `meta.stationxml`: Full StationXML record for this channel only.  Only
  present if keyword argument `full=true`.
"""
read(file, S::Type{<:Seis.Station{T}}; warn=true, full=false) where T =
    _parse_sxml(StationXML.read(file; warn), S; full)

"""
    parse(string, S::Type{Seis.Station}; warn=true, full=false) -> stations::Vector{S}

Parse a string containing a StationXML file, optionall
specifying a type `S` which must be a `Station` type.
By default, `S` is `GeogStation{Float64}`.
Return a vector of all channels specified in the StationXML file.

See [`SeisStationXML.read`](@ref) for more details.
"""
parse(string::AbstractString, S::Type{<:Seis.Station{T}}; warn=true, full=false) where T =
    _parse_sxml(StationXML.readstring(string; warn), S; full)

"""
    _parse_sxml(stationxml, S; warn=true, full=false) -> stations

See [`SeisStationXML.read`](@ref) for details.
"""
function _parse_sxml(
    stationxml::StationXML.FDSNStationXML,
    S::Type{<:Seis.Station{T}};
    warn=true,
    full=false
) where T
    stations = S[]
    for network in stationxml.network
        for station in network.station
            for channel in station.channel
                # Fill in information which may need adjusting first
                cha = S(
                    net=network.code, sta=station.code,
                    loc=channel.location_code, cha=channel.code,
                    # These things have a `.value` field containing the
                    # value of interest, so need special treatment
                    lon=_getifnotmissing(T, channel.longitude, :value),
                    lat=_getifnotmissing(T, channel.latitude, :value),
                    elev=_getifnotmissing(T, channel.elevation, :value),
                    azi=_getifnotmissing(T, channel.azimuth, :value),
                    inc=_getifnotmissing(T, channel.dip, :value)
                )

                # Then adjust information by converting to correct units, etc.,
                # which will preserve `missing` information
                cha.inc = cha.inc + 90
                cha.meta.burial_depth = _getifnotmissing(T, channel.depth, :value)
                cha.meta.startdate = coalesce(channel.start_date, station.start_date)
                cha.meta.enddate = coalesce(channel.end_date, station.end_date)
                if full
                    cha.meta.stationxml = filter_stationxml(
                        stationxml, network, station, channel
                    )
                end

                push!(stations, cha)
            end
        end
    end

    stations
end

"""
    filter_stationxml(sxml, network, station=nothing, channel=nothing) -> ::FDSNStationXML

Take a `FDSNStationXML` document and filter out all networks, stations
and channels apart from `network`, `station` and `channel`.
The returned copy only contains one network,
and maybe one station in that network if `station` is not `nothing`,
and maybe one channel at that station if `channel` is not `nothing`.
"""
function filter_stationxml(
    sxml::StationXML.FDSNStationXML, network, station=nothing, channel=nothing
)
    sta = _filtered_copy(station, channel)
    net = _filtered_copy(network, sta)
    _filtered_copy(sxml, net)
end

"""
    _filtered_copy(sta::StationXML.Station, cha::StationXML.Channel) -> ::StationXML.Station

Return a copy of `sta`, but with a single channel `cha` rather than all
channels.  If `sta.selected_num_channels` is defined, then the copy
will set this field to 1; otherwise it remains `missing`.
"""
function _filtered_copy(sta::StationXML.Station, cha::StationXML.Channel)
    sta_kwargs = Dict(
        field=>getproperty(sta, field) for field in propertynames(sta)
        if field !== :channel && field !== :selected_number_channels
    )
    num_chans = ismissing(sta.selected_number_channels) ? missing : 1
    StationXML.Station(; sta_kwargs..., channel=[cha], selected_number_channels=num_chans)
end

"""
    _filtered_copy(net::StationXML.Network, sta::StationXML.Station) -> ::StationXML.Network

Return a copy of `net`, but with a single station `sta` rather than all
stations.  If `sta.selected_num_stations` is defined, then the copy
will set this field to 1; otherwise it remains `missing`.
"""
function _filtered_copy(net::StationXML.Network, sta::StationXML.Station)
    net_kwargs = Dict(
        field=>getproperty(net, field) for field in propertynames(net)
        if field !== :station && field !== :selected_number_stations
    )
    num_stas = ismissing(net.selected_number_stations) ? missing : 1
    StationXML.Network(; net_kwargs..., station=[sta], selected_number_stations=num_stas)
end

"""
    _filtered_copy(sxml::StationXML.FDSNStationXML, net::StationXML.Network) -> ::StationXML.FDSNStationXML

Return a copy of `sxml`, but with a single network `net` rather than all
networks.
"""
function _filtered_copy(sxml::StationXML.FDSNStationXML, net::StationXML.Network)
    sxml_kwargs = Dict(
        field=>getproperty(sxml, field) for field in propertynames(sxml)
        if field !== :network
    )
    StationXML.FDSNStationXML(; sxml_kwargs..., network=[net])
end

"Return `missing` if `val is `missing`, and `val.field` otherwise."
_getifnotmissing(val, field) = val === missing ? missing : getfield(val, field)
"Return `missing` if `val` is missing, and `T(val.field)` otherwise."
_getifnotmissing(T, val, field) = val === missing ? missing : T(getfield(val, field))

end # module
