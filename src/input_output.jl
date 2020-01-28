"""
    read_sac(file) → ::Trace

Read a single evenly-sampled SAC file and return a Trace.

    read_sac(glob, dir) → ::Vector{Trace}

Read SAC files which match the patern `glob` in directory `dir` and return
a set of `Traces`.  Add the file names to `t.meta.file`.  These are relative
paths.

When reading SAC files, the following conventions are observed:

- The event id is held in KEVNM
- Channel ID is held in KCMPNM
- Location ID is held in KHOLE
- If O and the file origin time parameters are set, O is shifted to 0 time, and
  all time picks are adjusted.  This is similar to using the commands `ch o gmt [date];
  ch allt (0 - &1,o&)` to set the origin.
- Time picks are added to the `Trace` picks.

SAC headers which don't directly translate to `Trace` attributes are placed in the
.meta field and have names prefixed by "SAC_".
"""
function read_sac(glob, dir; kwargs...)
    s, f = SAC.read_wild(glob, dir; kwargs...)
    length(s) == 0 && return Trace{Float32,Vector{Float32},String}[]
    t = Trace.(s)
    t.meta.file = f
    t
end

function read_sac(file; kwargs...)
    t = Trace(SAC.read(file; kwargs...))
    t.meta.file = file
    t
end

"""
    Trace(s::SACtr) -> t

Construct the `Trace` `t` from the `SACtr` `s`.  See [read_sac](@ref) for details
of which headers are transferred to which fields in `t`.
"""
function Trace(s::SAC.SACtr)
    sac_trace_hdr = (:b, :e, :o, :npts, :delta, :depmin, :depmax, :depmen, :nvhdr, :leven)
    sac_evt_hdr = (:evlo, :evla, :evdp, :kevnm)
    sac_sta_hdr = (:stlo, :stla, :stel, :kstnm, :knetwk, :khole)
    sac_cmp_hdr = (:kcmpnm, :cmpaz, :cmpinc)
    sac_date_hdr = (:nzyear, :nzjday, :nzhour, :nzmin, :nzsec, :nzmsec)
    sac_time_hdr = (:a, :f, :b, :e, :t0, :t1, :t2, :t3, :t4, :t5, :t6, :t7, :t8, :t9)
    sac_picks_time_hdr = (:a, :f, :t0, :t1, :t2, :t3, :t4, :t5, :t6, :t7, :t8, :t9)
    sac_picks_name_hdr = (:ka, :kf, :kt0, :kt1, :kt2, :kt3, :kt4, :kt5, :kt6, :kt7, :kt8, :kt9)
    sac_derived_hdr = (:az, :baz, :gcarc, :dist)
    # All the headers we incorporate
    sac_hdr = (sac_trace_hdr..., sac_evt_hdr..., sac_sta_hdr..., sac_cmp_hdr...,
        sac_date_hdr..., sac_time_hdr..., sac_picks_time_hdr..., sac_picks_name_hdr...,
        sac_derived_hdr...)
    t = Trace{Float32}(s[:b], s[:delta], s[:t])
    # Origin time
    if !any(SAC.isundefined, getfield.(s, sac_date_hdr))
        t.evt.time = Dates.DateTime(s[:nzyear], 1, 1, s[:nzhour], s[:nzmin],
            s[:nzsec], s[:nzmsec]) + Dates.Day(s[:nzjday] - 1)
        # Seis.Traces are always aligned so that 0 is the origin, whatever it is
        # Shift all SAC headers such that O is 0.
        if !SAC.isundefined(s, :o)
            for h in sac_time_hdr
                !SAC.isundefined(s, h) && (s[h] -= s[:o])
            end
            t.b -= s[:o]
            s[:o] = 0
        end
    end

    # Event
    for (sacfield, tfield) in (:evlo=>:lon, :evla=>:lat, :evdp=>:dep, :kevnm=>:id)
        setproperty!(t.evt, tfield, _sacmissing(s[sacfield]))
    end

    # Station and component
    for (sacfield, tfield) in (:stlo=>:lon, :stla=>:lat)
        setproperty!(t.sta, tfield, _sacmissing(s[sacfield]))
    end
    s[:stel] !== missing && (t.sta.pos.dep = -s[:stel]/1000)
    for (sacfield, tfield) in (:kstnm=>:sta, :knetwk=>:net, :khole=>:loc,
                               :kcmpnm=>:cha, :cmpaz=>:azi, :cmpinc=>:inc)
        setproperty!(t.sta, tfield, _sacmissing(s[sacfield]))
    end

    # Time picks
    for (time_field, name_field) in zip(sac_picks_time_hdr, sac_picks_name_hdr)
        if !SAC.isundefined(s, time_field)
            key = Symbol(uppercase(String(time_field)))
            time = s[time_field]
            name = _sacmissing(s, name_field)
            t.picks[key] = Pick{Float32,String}((time, name))
        end
    end

    # Other headers
    for sacfield in SAC.sac_all_hdr
        sacfield in sac_hdr && continue
        SAC.isundefined(s, sacfield) && continue
        metafield = Symbol("SAC_" * String(sacfield))
        t.meta[metafield] = s[sacfield]
    end
    t
end

"""
    write_sac(t, file; littleendian=false)

Write the `Trace` `t` to `file` in SAC format.

Keys in the `t.meta` field which begin with `SAC_` have their values written to
the corresponding SAC field (e.g., `t.meta.SAC_kuser0` is written to the `KUSER0`
header).  The user is responsible for ensuring that the values corresponding to these
keys can be converted to the correct header type.  Note also that `SAC_` `meta` fields
override the equivalent `Trace` headers (e.g., `t.sta.sta` is equivalent to `SAC_kstnm`)
and so one way to override the values in `Trace` headers is to set the `SAC_` fields.

Time picks with keys corresponding to SAC picks headers (`A`, `F`, and `T0` to `T9`)
are transferred, but other picks are not.

By default, files are written to disk in bigendian format (MacSAC or SAC/BRIS
convention).  Use `littleendian=true` to write in littleendian byte order
(SAC/IRIS or SAC2000 convention).
"""
write_sac(t::AbstractTrace, file; littleendian=false) =
    SAC.write(SACtr(t), file; byteswap=!littleendian)

"""
    SACtr(t::Trace) -> s

Construct a `SACtr` s from s `Seis.Trace`.
"""
function SACtr(t::AbstractTrace)
    s = SACtr(trace(t), t.delta, t.b)
    for (sacfield, val) in (
                :o => 0,
                :stlo => t.sta.lon,
                :stla => t.sta.lat,
                :stel => t.sta.elev,
                :evlo => t.evt.lon,
                :evla => t.evt.lat,
                :evdp => t.evt.dep,
                :cmpaz => t.sta.azi,
                :cmpinc => t.sta.inc,
                :nzyear => ismissing(t.evt.time) ? nothing : Dates.year(t.evt.time),
                :nzjday => ismissing(t.evt.time) ? nothing : Dates.dayofyear(t.evt.time),
                :nzhour => ismissing(t.evt.time) ? nothing : Dates.hour(t.evt.time),
                :nzmin => ismissing(t.evt.time) ? nothing : Dates.minute(t.evt.time),
                :nzsec => ismissing(t.evt.time) ? nothing : Dates.second(t.evt.time),
                :nzmsec => ismissing(t.evt.time) ? nothing : Dates.millisecond(t.evt.time),
                :leven => true,
                :kstnm => t.sta.sta,
                :kevnm => t.evt.id,
                :khole => t.sta.loc,
                :kcmpnm => t.sta.cha,
                :knetwk => t.sta.net
            )
        !ismissing(val) && val !== nothing && (s[sacfield] = val)
    end
    # Time picks
    for (key, (time, name)) in t.picks
        keystring = string(key)
        # A, F or Tn marker
        if occursin(r"^([AF]|T[0-9])$"i, keystring)
            sac_tfield = Symbol(lowercase(keystring))
            s[sac_tfield] = time
            if !ismissing(name)
                sac_kfield = Symbol("k", sac_tfield)
                s[sac_kfield] = name
            end
        end
    end
    # Remaining headers
    for (name, val) in t.meta
        name_string = String(name)
        if occursin(r"^SAC_.*", name_string)
            sacfield = Symbol(name_string[5:end])
            s[sacfield] = val
        end
    end
    SAC.update_headers!(s)
    s
end

"""
    channel_code(t::Trace) -> code
    channel_code(s::Station) -> code

Return the channel code for trace `t` or station `s`, in the form of
"\$net.\$name.\$location.\$component".  Missing fields are left blank.
"""
channel_code(sta::Station) = join(_blankmissing.((sta.net, sta.sta, sta.loc, sta.cha)), ".")
channel_code(t::AbstractTrace) = channel_code(t.sta)

_blankmissing(x) = string(ismissing(x) ? "" : x)

_sacmissing(x) = SAC.isundefined(x) ? missing : x
_sacmissing(x::String) = (x = replace(x, "\0"=>""); SAC.isundefined(x) ? missing : x)
_sacmissing(s::SACtr, x) = _sacmissing(s[x])
