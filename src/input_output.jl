# TODO: Refactor all this so IO is its own module and formats are in separate
#       files.  Only export the functionality needed

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
    sac_sta_hdr = (:stlo, :stla, :stel, :kstnm, :knetwk)
    sac_cmp_hdr = (:kcmpnm, :cmpaz, :cmpinc)
    sac_date_hdr = (:nzyear, :nzjday, :nzhour, :nzmin, :nzsec, :nzmsec)
    sac_time_hdr = (:a, :f, :b, :e, :t0, :t1, :t2, :t3, :t4, :t5, :t6, :t7, :t8, :t9)
    sac_picks_time_hdr = (:a, :f, :t0, :t1, :t2, :t3, :t4, :t5, :t6, :t7, :t8, :t9)
    sac_picks_name_hdr = (:ka, :kf, :kt0, :kt1, :kt2, :kt3, :kt4, :kt5, :kt6, :kt7, :kt8, :kt9)
    sac_derived_hdr = (:az, :baz, :gcarc, :dist)
    # All the headers we incorporate
    sac_hdr = (sac_trace_hdr..., sac_evt_hdr..., sac_sta_hdr..., sac_cmp_hdr...,
        sac_date_hdr..., sac_time_hdr..., sac_picks_time_hdr..., sac_picks_time_hdr...,
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
        setfield!(t.evt, tfield, _sacmissing(s[sacfield]))
    end

    # Station and component
    for (sacfield, tfield) in (:stlo=>:lon, :stla=>:lat, :stel=>:elev,
                               :kstnm=>:sta, :knetwk=>:net, :khole=>:loc,
                               :kcmpnm=>:cha, :cmpaz=>:azi, :cmpinc=>:inc)
        setfield!(t.sta, tfield, _sacmissing(s[sacfield]))
    end

    # Time picks
    for (time_field, name_field) in zip(sac_picks_time_hdr, sac_picks_name_hdr)
        if !SAC.isundefined(s, time_field)
            name = if time_field in (:a, :f)
                SAC.isundefined(s, name_field) ?
                    uppercase(String(time_field)) : s[name_field]
            else
                _sacmissing(s, name_field)
            end
            add_pick!(t, s[time_field], name)
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
    write_sac(t, file)

Write the `Trace` `t` to `file` in SAC big-endian format.

Keys in the `t.meta` field which begin with "SAC_" have their values written to
the corresponding SAC field (e.g., `t.meta[:SAC_kuser0]` is written to the KUSER0
header).  The user is responsible for ensuring that the values corresponding to these
keys can be converted to the correct header type.
"""
write_sac(t::AbstractTrace, file) = SAC.write(SACtr(t), file)

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
    ipick = 0
    for (time, name) in picks(t)
        if occursin(r"^[AF]$"i, name)
            sac_tfield = Symbol(lowercase(name))
            s[sac_tfield] = time
        else
            if ipick > 9
                @warn("Only the first 10 picks added to the `SACtr`")
                break
            end
            sac_tfield = Symbol(:t, ipick)
            sac_kfield = Symbol(:kt, ipick)
            s[sac_tfield] = time
            s[sac_kfield] = name
            ipick += 1
        end
    end
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
    parse_mseed([T=$(DEFAULT_FLOAT)[[, V=Vector{T},] S=$DEFAULT_STRING]], data) -> t::Vector{Trace{T,V,S}}

Parse the contents of `data` (usually a `Vector{UInt8}`) as miniSEED-formatted
seismic data, and return the traces `t`.

Optionally specify the parameters of the `Trace{T,V,S}` returned with `T`, `V` and `S`.

miniSEED-specific headers are stored in `t.meta` with names beginning `MSEED_`.
"""
parse_mseed(T, V, S, data) = (t=read_mseed(T, V, S, IOBuffer(data)); t.meta.file=missing; t)
parse_mseed(T, S, data) = parse_mseed(T, Vector{T}, S, data)
parse_mseed(T, data) = parse_mseed(T, Vector{T}, DEFAULT_STRING, data)
parse_mseed(data) = parse_mseed(DEFAULT_FLOAT, Vector{DEFAULT_FLOAT}, String, data)

"""
    read_mseed([T=$DEFAULT_FLOAT[, V=Vector{T}, S=$DEFAULT_STRING]], file) -> t::Vector{Trace{T,V,S}}

Read a single miniSEED-formatted `file` and return a set of `Trace`s `t` containing
the contents.

Optionally specify the parameters of the `Trace{T,V,S}` returned with `T`, `V` and `S`.

miniSEED-specific headers are stored in `t.meta` with names beginning `MSEED_`.
"""
function read_mseed(T, V, S, file)
    stream = ObsPy[:read](file)
    t = Trace{T,V,S}[]
    for (i, s) in enumerate(stream)
        b = s[:times]()[1]
        delta = s[:stats][:delta]
        trace = convert(V, s[:data])
        tt = Trace{T,V,S}(b, delta, trace)
        tt.sta.net = s[:stats][:network]
        tt.sta.sta = s[:stats][:station]
        tt.sta.loc = s[:stats][:location]
        tt.sta.cha = s[:stats][:channel]
        # FIXME: Use UTCDateTime epoch of 1970-01-01T00:00:00.000000000
        #        Not clear how to use public interface via PyCall
        tt.evt.time = DateTime(1970) +
            Dates.Nanosecond(s[:stats][:starttime][:_UTCDateTime__ns])
        # Add miniSEED-specific flags
        # FIXME: Is there a better way in Julia to get the keys in a Python dict?
        for (k, v) in s[:stats][:mseed][:__dict__]
            k′ = Symbol(string("MSEED_", k))
            tt.meta[k′] = v
        end
        tt.meta.file = file
        push!(t, tt)
    end
   t
end

read_mseed(T, S, file) = read_mseed(T, Vector{T}, S, file)
read_mseed(T, file) = read_mseed(T, Vector{T}, DEFAULT_STRING, file)
read_mseed(file) = read_mseed(DEFAULT_FLOAT, Vector{DEFAULT_FLOAT}, DEFAULT_STRING, file)

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
_sacmissing(s::SACtr, x) = _sacmissing(s[x])
