# Event-station geometry calculations on a (maybe flattened) sphere

"""
    azimuth(trace; sphere=false, flattening=Geodesics.F_WGS84) -> az
    azimuth(event, station; sphere=false, flattening=Geodesics.F_WGS84) -> az

Return the azimuth `az` from the event to the station (a seismic station) in
degrees east from local north at the event for a `trace`.  Alternatively specify
the `event` and `station` individually.

Optionally specify the `flattening` of the ellipsoid of rotation on which this
is computed, which defaults to that of the WGS84 ellipsoid.  If `sphere` is `true`,
then `flattening` is set to zero and the calculation is performed on a sphere.
"""
function azimuth(e::GeogEvent, s::GeogStation; sphere=false, flattening=Geodesics.F_WGS84)
    _check_headers_geometry(e, s)
    Geodesics.azimuth(e.lon, e.lat, s.lon, s.lat, true, f=sphere ? 0.0 : flattening)
end

function azimuth(e::CartEvent, s::CartStation)
    any(x->x===missing, (e.x, e.y, s.x, s.y)) &&
        throw(ArgumentError("Insufficient information in trace to compute geometry." *
                            " (Need: evt.x, evt.y, sta.x and sta.y"))
    mod(rad2deg(atan(s.x - e.x, s.y - e.y)), 360)
end

azimuth(t::AbstractTrace, args...; kwargs...) = azimuth(t.evt, t.sta, args...; kwargs...)

"""
    backazimuth(trace; flattening=$(Geodesics.F_WGS84)) -> baz
    backazimuth(station, event; flattening=$(Geodesics.F_WGS84)) -> baz

Return the backazimuth `baz` from the station (a seismic station) to an event in
degrees east from local north at the station for a `trace`.  Alternatively specify
the `station` and `event` individually.

Optionally specify the `flattening` of the ellipsoid of rotation on which this
is computed, which defaults to that of the WGS84 ellipsoid.  If `sphere` is `true`,
then `flattening` is set to zero and the calculation is performed on a sphere.
"""
function backazimuth(s::GeogStation, e::GeogEvent; sphere=false,
        flattening=Geodesics.F_WGS84)
    _check_headers_geometry(e, s)
    Geodesics.azimuth(s.lon, s.lat, e.lon, e.lat, true, f=(sphere ? 0.0 : flattening))
end
backazimuth(t::AbstractTrace, args...; kwargs...) = backazimuth(t.sta, t.evt; kwargs...)

backazimuth(s::CartStation, e::CartEvent) = mod(azimuth(e, s) + 180, 360)

"""
    incidence(event, station) -> i
    incidence(trace) -> i

Return the angle of incidence `i`° between a cartesian `event` and `station`,
or a cartesian `trace`.  The angle of incidence is defined downwards from the
positive z (upward) direction.
"""
function incidence(e::CartEvent, s::CartStation)
    h = distance_km(e, s)
    rad2deg(atan(h, (s.z - e.z)/1000))
end
incidence(t::AbstractTrace) = incidence(t.evt, t.sta)

"""
    distance_deg(trace; sphere=false, flattening=$(Geodesics.F_WGS84)) -> Δ
    distance_deg(station, event; sphere=false, flattening=$(Geodesics.F_WGS84)) -> Δ
    distance_deg(event, station; sphere=false, flattening=$(Geodesics.F_WGS84)) -> Δ

For a `trace` or an `event`-`station` pair, return the epicentral angular distance
`Δ` in degrees.

Optionally specify the `flattening` of the ellipsoid of rotation on which this
is computed, which defaults to that of the WGS84 ellipsoid.  If `sphere` is `true`,
then `flattening` is set to zero and the calculation is performed on a sphere.
"""
function distance_deg(e::GeogEvent, s::GeogStation;
                      sphere=false, flattening=Geodesics.F_WGS84)
    _check_headers_geometry(e, s)
    Geodesics.angular_distance(e.lon, e.lat, s.lon, s.lat, true, f=(sphere ? 0.0 : flattening))
end
distance_deg(s::Station, e::Event; kwargs...) = distance_deg(e, s; kwargs...)
distance_deg(t::Trace; kwargs...) = distance_deg(t.evt, t.sta; kwargs...)

"""
    distance_km(event, station; sphere=false, a=$(Geodesics.EARTH_R_MAJOR_WGS84/1e3), flattening=$(Geodesics.F_WGS84)) -> d

For a geographic `trace` or `event`–`station` pair, return the epicentral
surface distance `d` in km between them.

Optionally specify the semimajor radius `a` in km and `flattening` of the
ellipsoid of rotation on which this is computed, which defaults to that of the
WGS84 ellipsoid.    If `sphere` is `true`, then `flattening` is set to zero and
the calculation is performed on a sphere.
"""
function distance_km(e::GeogEvent, s::GeogStation;
                     sphere=false, a=Geodesics.EARTH_R_MAJOR_WGS84/1e3,
                     flattening=Geodesics.F_WGS84)
    _check_headers_geometry(e, s)
    Geodesics.surface_distance(e.lon, e.lat, s.lon, s.lat, a, true,
                               f=(sphere ? 0.0 : flattening))
end

"""
    distance_km(event, station) -> d

For a cartesian `trace` or `event`-`station` pair, return the
epicentral distance in km between them.
"""
function distance_km(e::CartEvent, s::CartStation)
    _check_headers_geometry(e, s)
    any(ismissing, (e.z, s.z)) &&
        throw(ArgumentError("z coordinates must be defined for `Seis.Cartesian` points"))
    sqrt((e.x - s.x)^2 + (e.y - s.y)^2)/1000
end

distance_km(s::Station, e::Event; kwargs...) = distance_km(e, s; kwargs...)
distance_km(t::AbstractTrace; kwargs...) = distance_km(t.evt, t.sta; kwargs...)

"""
    _check_headers_geometry(evt, sta) -> nothing

Check that an `evt`–`sta` pair have all the necessary information to compute
the azimuth, backazimuth and epicentral distance.  If they do not, throw an
`ArgumentError`.
"""
function _check_headers_geometry(evt::GeogEvent, sta::GeogStation)
    any(ismissing, (evt.lon, evt.lat, sta.lon, sta.lat)) &&
        throw(ArgumentError("Insufficient information in trace to compute geometry." *
                            " (Need: evt.lon, evt.lat, sta.lon, sta.lat)")) ||
    nothing
end
function _check_headers_geometry(evt::CartEvent, sta::CartStation)
    any(ismissing, (evt.x, evt.y, evt.z, sta.x, sta.y, sta.z)) &&
        throw(ArgumentError("Insufficient information in trace to compute geometry." *
                            " (Need: evt.x, evt.y, evt.z, sta.x, sta.y, sta.z"))
    nothing
end
