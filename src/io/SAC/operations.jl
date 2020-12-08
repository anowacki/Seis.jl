"""
    update_headers!(s::SACTrace)

Ensure that header values which are based on the trace or other header values
are consistent, such as `depmax`.  Should be called after any operation on the trace
`s.t`.

!!! note
    If the trace is empty, then `depmax`, `depmin` and `depmen` are not updated,
    but the trace end time is adjusted according to `npts`.
"""
function update_headers!(s::SACTrace)
    if !isempty(s.t)
        s.depmax = maximum(s.t)
        s.depmin = minimum(s.t)
        s.depmen = mean(s.t)
    end
    s.e = s.b + s.delta*(s.npts - 1)
    s
end

"""
    update_great_circle!(s::SACtr)

If all headers `evlo`, `evla`, `stlo` and `stla` are set, update the values of
`az`, `baz` and `gcarc`.
"""
function update_great_circle!(s::SACTrace)
    any([s.evlo, s.evla, s.stlo, s.stla] .== SAC_RNULL) && return
    s.az = Geodesics.azimuth(s.evlo, s.evla, s.stlo, s.stla, true,
                             f=Geodesics.F_WGS84)
    s.baz = Geodesics.azimuth(s.stlo, s.stla, s.evlo, s.evla, true,
                              f=Geodesics.F_WGS84)
    s.gcarc = Geodesics.angular_distance(s.evlo, s.evla, s.stlo, s.stla, true,
                                         f=Geodesics.F_WGS84)
    s
end
