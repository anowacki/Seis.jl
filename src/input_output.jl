#
# SAC
#

"""
    read_sac(file; terse=false, header_only=true) → ::Trace

Read a single evenly-sampled SAC file and return a Trace.  If `terse` is
`true`, then warn when auto-byteswapping files.  To read only SAC
headers from files, returning an empty trace, set `header_only` to `true.

---

    read_sac(glob, dir; echo=true, header_only=false) → ::Vector{Trace}

Read SAC files which match the patern `glob` in directory `dir` and return
a set of `Traces`.  Add the file names to `t.meta.file`.  These are relative
paths.

File names matching the pattern are shown unless `echo` is `false`.

---

When reading SAC files, the following conventions are observed:

- The event id is held in header `KEVNM`
- Channel ID is held in `KCMPNM`
- Location ID is held in `KHOLE`
- If `O` and the file origin time parameters are set, `O` is shifted to 0 time, and
  all time picks are adjusted.  This is similar to using the commands `ch o gmt [date];
  ch allt (0 - &1,o&)` to set the origin in SAC.
- Time picks are added to the `Trace` picks.

SAC headers which don't directly translate to `Trace` attributes are placed in the
`.meta` field and have names prefixed by `"SAC_"`.
"""
function read_sac(glob, dir; kwargs...)
    s, f = SAC.read_wild(glob, dir; kwargs...)
    length(s) == 0 && return Trace{Float32, Vector{Float32}, Geographic{Float32}}[]
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
    Trace(s::SAC.SACTrace) -> t

Construct the `Trace` `t` from the `SACtr` `s`.  See [read_sac](@ref) for details
of which headers are transferred to which fields in `t`.
"""
function Trace(s::SAC.SACTrace)
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
    !SAC.isundefined(s, :stel) && (t.sta.pos.dep = -s[:stel]/1000)
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
            t.picks[key] = Pick{Float32}(time, name)
        end
    end

    # Other headers
    for sacfield in SAC.SAC_ALL_HDR
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
Note that the header is lowercase (i.e., `SAC_kstnm` not `SAC_KSTNM`).

Time picks with keys corresponding to SAC picks headers (`A`, `F`, and `T0` to `T9`)
are transferred, but other picks are not.

If `t` is in a Cartesian reference frame (i.e., its positions are given
by `CartEvent` and `CartStation`), then the Cartesian station coordinates
`x`, `y` and `z` are saved respectively to headers `USER0`, `USER1` and `USER2`.
Likewise, the event coordinates are saved respectively to `USER3`, `USER4` and
`USER5`.  Any information in `meta` fields `SAC_user0` to `SAC_user5` will
overwrite this data.

!!! note
    The convention on how non-geographic coordinates are written in SAC headers is
    not part of the API and may change at any time.  Saving non-standard information
    in SAC headers should be done explicitly by the user if this information is
    important.

By default, files are written to disk in bigendian format (MacSAC or SAC/BRIS
convention).  Use `littleendian=true` to write in littleendian byte order
(SAC/IRIS or SAC2000 convention).
"""
write_sac(t::AbstractTrace, file; littleendian=false) =
    SAC.write(SAC.SACTrace(t), file; byteswap=!littleendian)

"""
    SACTrace(t::Trace) -> s

Construct a `SACTrace` `t` from a `Trace` `s`.
"""
function SAC.SACTrace(t::AbstractTrace)
    s = SAC.SACTrace(trace(t), t.delta, t.b)
    position_headers = _sac_position_headers(t)
    for (sacfield, val) in pairs((
                o = 0,
                position_headers...,
                cmpaz = t.sta.azi,
                cmpinc = t.sta.inc,
                nzyear = ismissing(t.evt.time) ? nothing : Dates.year(t.evt.time),
                nzjday = ismissing(t.evt.time) ? nothing : Dates.dayofyear(t.evt.time),
                nzhour = ismissing(t.evt.time) ? nothing : Dates.hour(t.evt.time),
                nzmin = ismissing(t.evt.time) ? nothing : Dates.minute(t.evt.time),
                nzsec = ismissing(t.evt.time) ? nothing : Dates.second(t.evt.time),
                nzmsec = ismissing(t.evt.time) ? nothing : Dates.millisecond(t.evt.time),
                leven = true,
                kstnm = t.sta.sta,
                kevnm = t.evt.id,
                khole = t.sta.loc,
                kcmpnm = t.sta.cha,
                knetwk = t.sta.net,
                iftype = SAC.SAC_ITIME
            ))
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

"Fields relating to positions which are set differently depending
on the geometry of the traces."
_sac_position_headers(t::Trace{T,V,P}) where {T,V,P<:Geographic}= (
    stlo = t.sta.lon,
    stla = t.sta.lat,
    stel = t.sta.elev,
    evlo = t.evt.lon,
    evla = t.evt.lat,
    evdp = t.evt.dep
    )
function _sac_position_headers(t::Trace{T,V,P}) where {T,V,P<:Cartesian}
    any_present = any(x -> x !== missing, (t.sta.x, t.sta.y, t.sta.z, t.evt.x, t.evt.y, t.evt.z))
    (
    user0 = t.sta.x,
    kuser0 = any_present ? "Seis.jl" : nothing,
    user1 = t.sta.y,
    kuser1 = any_present ? "xyz pos" : nothing,
    user2 = t.sta.z,
    user3 = t.evt.x,
    user4 = t.evt.y,
    user5 = t.evt.z,
    )
end
_sac_position_headers(_::AbstractTrace) = ()

"Convert a SAC undefined header value into `missing` or return its value"
_sacmissing(x) = SAC.isundefined(x) ? missing : x
# Remove trailing, illegal null characters from strings; the spec is to pad with spaces.
_sacmissing(x::String) = (x = replace(x, "\0"=>""); SAC.isundefined(x) ? missing : x)
_sacmissing(s::SAC.SACTrace, x) = _sacmissing(s[x])

#
# Miniseed
#

"""
    read_mseed(file; kwargs...) -> traces
    read_mseed(file, T; kwargs...) -> traces::Vector{T}

Read a single miniseed file from disk and return a set of `Trace`s.

The `meta.mseed_file` field of each trace contains the file name. 

Optionally specify the type of trace `T <: AbstractTrace` to read.  By default,
`T` is `$(Miniseed.DEFAULT_TRACE)`, since almost all
seismic data stored in Miniseed format is single-precision, and because
the sampling rate is stored at a 64-bit float in miniSEED files.

# Example
Read a single file:
```
julia> read_mseed("data.mseed")
```

Read a single file assuming a Cartesian geometry:
```
julia> read_mseed(CartTrace{Float64, Vector{Float32}}, "data.mseed")
```


# Handling gapped/overlapped data

When channels containing gaps or overlaps are encountered, they are
split into multiple `Trace`s as each `Trace` must be continuous and evenly
sampled.  However, data quite often contain single-sample offsets which
are later corrected, and so these are ignored by default.

Use the keyword argument `maximum_gap` to control
whether or not gaps cause new traces to be created.  See below for more details.


# Keyword arguments

The following keyword arguments can be passed to `read_mseed`:

- `maximum_gap`: The maximum absolute gap length in s beyond which gaps are
  no longer tolerated in a single trace.  By default this is the sampling
  interval of the trace being read.

  !!! note
      Set `maximum_gap` to 0 to always split miniseed files into separate
      traces at all gaps.

- `verbose`: An integer starting from 0 upwards indicating how much
  information about the reding process should be printed to the screen.
  The default (`0`) only produces output for errors and warnings.

---

    read_mseed(data::Vector{UInt8}[, T]; kwargs...) -> traces

Read Miniseed `data` from memory, held as a set of bytes, optionally specifying
the type `T` of traces to return.  Keyword arguments are the same as for
reading from a file on disk.
"""
read_mseed(file, T::Type=Miniseed.DEFAULT_TRACE; kwargs...) =
    Miniseed.read(file, T; kwargs...)

"""
    read_mseed(pattern, dir) -> ::Vector{<:Trace}
    read_mseed(pattern, dir, T) -> Vector{T}

Read all files matching `pattern` in directory `dir`.

See `Glob.glob` for details of pattern matching.

Optionally specify the type of trace `T <: AbstractTrace` to read.

# Example
Read all files matching `"TA.*.BHZ.mseed"` in all directories within
`DATA` which themselves match `"Event_??"`:
```
julia> read_mseed("Event_??/TA.*.BHZ.mseed", "DATA")
```
"""
function read_mseed(pattern, dir, T::Type=Miniseed.DEFAULT_TRACE; kwargs...)
    files = Glob.glob(pattern, dir)
    if isempty(files)
        T[]
    else
        reduce(vcat, [read_mseed(file, T; kwargs...) for file in files])
    end
end

"""
    write_mseed(file, t; append=false, verbose=0, pubversion=1, record_length=nothing, version=2)

Write the data contained in the`t` to `file` on disk
in miniSEED format.

`t` may be either a single `AbstractTrace`, or an array of `AbstractTrace`s.

miniSEED files can contain data in (amongst others) in `Float32`, `Float64`
and `Int32` format.  This function will use whatever precision or type the
data in `t` have and attempt to write.  If for example you want to write
a trace with a `Float64` element type to a miniSEED file with element
type `Float32`, you should first convert the trace using `convert`.
(Note that since Seis does not support integer-valued trace data, it will
not write 32-bit integer miniSEED files.)

If the trace does not have an origin time set, an error is thrown.

# Keyword arguments

- `append::Bool`: If `true`, add the data in `t` to the end of any data already
  existing in `file`.  miniSEED files can contain multiples traces.
- `verbose::Integer`: Controls the verbosity of the miniSEED conversion and writing
  process.  Larger values of `verbose` cause more output to be produced.
- `pubversion::Integer`: The publication version of data describes whether a
  set of data has been updated since being initially published.  Higher numbers
  correspond to records which supercede lower versions, which start at 1 (the default).
- `record_length::Integer`: The number of bytes used to write each miniSEED record.
- `version`: The miniSEED file version to write.  Can be `2` (the default) or `3`.
"""
write_mseed(file, t; kwargs...) = Miniseed.write(file, t; kwargs...)

"""
    channel_code(t::Trace) -> code
    channel_code(s::Station) -> code

Return the channel code for trace `t` or station `s`, in the form of
`"⟨network⟩.⟨name⟩.⟨location⟩.⟨component⟩"`.  Missing fields are left blank.
The information is taken respectively from the `net`, `sta`, `cha` and `loc`
fields of the `Station`.
"""
channel_code(sta::Station) = join(_blankmissing.((sta.net, sta.sta, sta.loc, sta.cha)), ".")
channel_code(t::AbstractTrace) = channel_code(t.sta)

_blankmissing(x) = string(ismissing(x) ? "" : x)
