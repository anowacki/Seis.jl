"""
    azimuth(trace; flattening=Geodesics.F_WGS84) -> az
    azimuth(event, station; flattening=Geodesics.F_WGS84) -> az

Return the azimuth `az` from the event to the station (a seismic station) in
degrees east from local north at the event for a `trace`.  Alternatively specify
the `event` and `station` individually.

Optionally specify the `flattening` of the ellipsoid of rotation on which this
is computed, which defaults to that of the WGS84 ellipsoid.
"""
azimuth(e::Event, s::Station; flattening=Geodesics.F_WGS84) =
    Geodesics.azimuth(e.lon, e.lat, s.lon, s.lat, true, f=flattening)
azimuth(t::Trace, args...) = azimuth(t.evt, t.cha, args...)

"""
    backazimuth(trace; flattening=$(Geodesics.F_WGS84)) -> baz
    backazimuth(station, event; flattening=$(Geodesics.F_WGS84)) -> baz

Return the backazimuth `baz` from the station (a seismic station) to an event in
degrees east from local north at the station for a `trace`.  Alternatively specify
the `station` and `event` individually.

Optionally specify the `flattening` of the ellipsoid of rotation on which this
is computed, which defaults to that of the WGS84 ellipsoid.
"""
backazimuth(s::Station, e::Event; flattening=Geodesics.F_WGS84) =
    Geodesics.azimuth(s.lon, s.lat, e.lon, e.lat, true, f=flattening)
backazimuth(t::Trace, args...; kwargs...) = backazimuth(t.cha, t.evt, args...; kwargs...)

"""
    distance_deg(trace; flattening=$(Geodesics.F_WGS84)) -> Δ
    distance_deg(station, event; flattening=$(Geodesics.F_WGS84)) -> Δ
    distance_deg(event, station; flattening=$(Geodesics.F_WGS84)) -> Δ

For a `trace` or an `event`-`station` pair, return the epicentral angular distance
`Δ` in degrees.

Optionally specify the `flattening` of the ellipsoid of rotation on which this
is computed, which defaults to that of the WGS84 ellipsoid.
"""
distance_deg(e::Event, s::Station; flattening=Geodesics.F_WGS84) =
    Geodesics.angular_distance(e.lon, e.lat, s.lon, s.lat, true, f=flattening)
distance_deg(s::Station, e::Event, args...; kwargs...) = distance_deg(e, c, args...; kwargs...)
distance_deg(t::Trace, args...; kwargs...) = distance_deg(t.evt, t.sta, args...; kwargs...)

"""
    distance_km(event, station; a=$(Geodesics.EARTH_R_MAJOR_WGS84/1e3), flattening=$(Geodesics.F_WGS84)) -> d

For a `trace` or `event`–`station` pair, retun the epicentral surface distance `d`
in km.

Optionally specify the semimajor radius `a` in km and `flattening` of the
ellipsoid of rotation on which this is computed, which defaults to that of the
WGS84 ellipsoid.
"""
distance_km(e::Event, s::Station; a=Geodesics.EARTH_R_MAJOR_WGS84/1e3,
              flattening=Geodesics.F_WGS84) =
    Geodesics.surface_distance(e.lon, e.lat, s.lon, s.lat, a, true, f=flattening)
distance_km(s::Station, e::Event; kwargs...) = distance_km(e, s; kwargs...)
distance_km(t::Trace; kwargs...) = distance_km(t.evt, t.sta; kwargs...)
