"""
# `SeisStationXML`

Module for reading and writing StationXML files into Seis structures
"""
module SeisStationXML

import Dates
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
                cha.meta.startdate = station.start_date
                cha.meta.enddate = station.end_date
                cha.meta.channel_startdate = channel.start_date
                cha.meta.channel_enddate = channel.end_date
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
    write(file, stations::AbstractArray{<:Seis.GeogStation}; use_stationxml=true, kwargs...)
    write(io, stations::AbstractArray{<:Seis.GeogStation}; use_stationxml=true, kwargs...)

Save the station metadata available in `stations` as a single StationXML
file on disk as `file`, which may be a string or `IO` stream object.

See [`Seis.write_stationxml`](@ref) for details.
"""
function write(file, stations; use_stationxml=true, kwargs...)
    sxml = StationXML.FDSNStationXML(stations; use_stationxml, kwargs...)
    Base.write(file, sxml)
end

"""
    StationXML.FDSNStationXML(stations::AbstractArray{<:Seis.GeogStation}; use_stationxml=true, kwargs...)

Create a new instance of [`StationXML.FDSNStationXML`](@ref) from a set
of `Seis.GeogStation`s.

If `use_stationxml` is `true` (the default), take information from the
StationXML stored in the `.meta.stationxml` field of each station.  Otherwise,
use only the data held in the standard [`Seis.Station`](@ref) structure.

`kwargs` are passed to the standard `StationXML.FDSNStationXML` constructor
and therefore can be used to override the default values for attributes and
elements in the output.

!!! note
    Only `Seis.GeogStation`s are supported because StationXML requires
    that coordinates be given in longitude and latitude.
"""
function StationXML.FDSNStationXML(
    stations::AbstractArray{<:Seis.GeogStation};
    use_stationxml=true,
    kwargs...
)
    if isempty(stations)
        throw(ArgumentError(
            "cannot construct FDSNStationXML from an empty set of stations"
        ))
    end

    _check_station_missing_data_and_throw(stations)

    # Count the number of channels per station and stations per network
    stations_per_network = Dict{String,Set{String}}()
    channels_per_station = Dict{Tuple{String,String},Set{Tuple{String,String}}}()
    for station in stations
        @show Seis.channel_code(station)
        net, sta, loc, cha = station.net, station.sta, coalesce(station.loc, ""), station.cha

        if !haskey(stations_per_network, net)
            stations_per_network[net] = Set{String}()
        end
        push!(stations_per_network[net], sta)

        if !haskey(channels_per_station, (net, sta))
            channels_per_station[net,sta] = Set{Tuple{String,String}}()
        end
        push!(channels_per_station[net,sta], (loc, cha))
    end

    channels_sxml = [StationXML.Channel(sta; use_stationxml) for sta in stations]
    stations_sxml = [
        StationXML.Station(
            sta; use_stationxml, channel=[channel_sxml],
            selected_number_channels=length(channels_per_station[sta.net,sta.sta])
        )
        for (sta, channel_sxml) in zip(stations, channels_sxml)
    ]
    networks_sxml = [
        StationXML.Network(
            sta; use_stationxml, station=[station_sxml],
            selected_number_stations=length(stations_per_network[sta.net])
        )
        for (sta, station_sxml) in  zip(stations, stations_sxml)
    ]
    sxmls = [
        StationXML.FDSNStationXML(;
            source="Seis.jl",
            module_name="Seis.jl: write_stationxml",
            module_uri="https://github.com/anowacki/Seis.jl",
            created=Dates.now(),
            schema_version="1.1",
            network=[network_sxml],
            kwargs...
        )
        for (sta, network_sxml) in zip(stations, networks_sxml)
    ]

    sxml1, sxml_rest = Iterators.peel(sxmls)
    foldl(merge!, sxml_rest; init=sxml1)
end
StationXML.FDSNStationXML(stations::Union{Seis.CartStation,AbstractArray{<:Seis.CartStation}}; kwargs...) =
    throw(ArgumentError("cannot construct StationXML from stations with Cartesian coordinates"))
StationXML.FDSNStationXML(station::Seis.GeogStation; kwargs...) =
    StationXML.FDSNStationXML([station]; kwargs...)

"""
    _check_station_missing_data_and_throw(stations)

Throw an `ArgumentError` if any of `stations` does not contain sufficient
information to construct an [`StationXML.FDSNStationXML`](@ref) object
and report which field is missing.
"""
function _check_station_missing_data_and_throw(stations)
    for sta in stations
        for field in (:net, :sta, :cha, :azi, :inc, :lon, :lat)
            if ismissing(getproperty(sta, field))
                throw(ArgumentError(_station_missing_data_string(sta, field)))
            end
        end
    end
end
_station_missing_data_string(sta, field) =
    "station $(sta.net).$(sta.sta).$(sta.loc).$(sta.cha) is missing field '$field' required to form StationXML"

"""
    StationXML.Network(sta::Seis.Station; use_stationxml=true, kwargs...)

Create a [`StationXML.Network`](@ref) from a `Station`.  `kwargs` are passed
on to the `Network` constructor so can be used to override any values within
`sta` or `sta.meta.stationxml`,
"""
function StationXML.Network(sta::Seis.Station; use_stationxml=true, kwargs...)
    if use_stationxml && haskey(sta.meta, :stationxml)
        net = sta.meta.stationxml.network[1]
        StationXML.Network(;
            (
                field => getproperty(net, field)
                for field in propertynames(net)
            )...,
            code=sta.net,
            kwargs...
        )
    else
        StationXML.Network(;
            code=sta.net,
            kwargs...
        )
    end
end

"""
    StationXML.Station(sta::Seis.Station; use_stationxml=true, kwargs...)

Create a [`StationXML.Station`](@ref) from a `Seis.Station`.  `kwargs` are passed
on to the `StationXML.Station` constructor so can be used to override any
values within `sta` or `sta.meta.stationxml`,
"""
function StationXML.Station(sta::Seis.Station; use_stationxml=true, kwargs...)
    # Arguments which may be overridden by any StationXML info in
    # `.meta.stationxml` but which are necessary to make the struct.
    default_kwargs = (
        site = StationXML.Site(name=""),
    )
    # Arguments which overwrite the StationXML info
    sta_kwargs = (
        code=sta.sta,
        latitude=sta.lat,
        longitude=sta.lon,
        elevation=coalesce(sta.elev, 0),
        start_date=sta.meta.startdate,
        end_date=sta.meta.enddate,
    )

    if use_stationxml && haskey(sta.meta, :stationxml)
        sta_stationxml = sta.meta.stationxml.network[1].station[1]
        StationXML.Station(;
            default_kwargs...,
            (
                field => getproperty(sta_stationxml, field)
                for field in propertynames(sta_stationxml)
            )...,
            sta_kwargs...,
            kwargs...
        )
    else
        StationXML.Station(;
            default_kwargs...,
            sta_kwargs...,
            kwargs...
        )
    end
end

"""
    StationXML.Channel(sta::Seis.Station; use_stationxml=true, kwargs...)

Create a [`StationXML.Channel`](@ref) from a `Station`.  `kwargs` are passed
on to the `Channel` constructor so can be used to override any values within
`sta` or `sta.meta.stationxml`.

Note that `sta.meta.startdate` will be used if no channel start date is present,
and likewise for station and channel end dates.
"""
function StationXML.Channel(sta::Seis.Station; use_stationxml=true, kwargs...)
    cha_kwargs = (
        code = sta.cha,
        start_date = coalesce(sta.meta.channel_startdate, sta.meta.startdate),
        end_date = coalesce(sta.meta.channel_enddate, sta.meta.enddate),
        latitude = sta.lat,
        longitude = sta.lon,
        elevation = coalesce(sta.elev, 0),
        depth = coalesce(get(sta.meta, :burial_depth, 0), 0),
        azimuth = sta.azi,
        dip = sta.inc - 90,
        location_code = coalesce(sta.loc, ""),
    )

    if use_stationxml && haskey(sta.meta, :stationxml)
        cha = sta.meta.stationxml.network[1].station[1].channel[1]
        StationXML.Channel(;
            (
                field => getproperty(cha, field)
                for field in propertynames(cha)
            )...,
            cha_kwargs...,
            kwargs...
        )
    else
        StationXML.Channel(;
            cha_kwargs...,
            kwargs...
        )
    end
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
