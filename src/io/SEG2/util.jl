"""
    _parse_warn(T, dict, key, default=zero(T); trace_number=nothing, warn=true) where T -> value

Parse a string as `T` from the trace or file descriptor headers.
Pass in the `dict` whose `key` if present will be parsed as `T`, and
otherwise return `default`.
"""
function _parse_warn(
    ::Type{T}, dict, key, default=zero(T);
    trace_number=nothing, warn=true,
) where T
    if haskey(dict, key)
        value = tryparse(T, dict[key])
        if isnothing(value)
            warn && @warn "Cannot parse \"$(dict[key])\" as a $T for header \"$key\" for trace number $trace_number"
            default
        else
            value
        end
    else
        warn && @warn "No header \"$key\" for trace number $trace_number; using default value"
        default
    end
end

"""
    _parse_warn_location(T, meta::Seis.SeisDict, key::Symbol, unit::AbstractString; warn::Bool=true) where T -> x, y, z
    _parse_warn_location(T, string, unit; warn::Bool=true, key=Symbol("")) where T -> x, y, z

Try to get the position fields, where `key` is either `:receiver_position`
or `:source_position`.  Return three values, of which all may be `missing` if
the information is not present, only the first is not `missing` if only
the x coordinate is given, or all three may be non-`missing` if all
three coordinates are given.

Values are converted to m according to the `unit` string obtained from
the file descriptor.
"""
function _parse_warn_location(
    ::Type{T}, str::AbstractString, unit::AbstractString;
    key=:unspecified,
    trace_number=nothing,
    warn::Bool=true
) where T
    isempty(str) && return missing, missing, missing

    tokens = split(str)
    if length(tokens) >= 1
        x = tryparse(T, tokens[1])

        if length(tokens) >= 2
            y = tryparse(T, tokens[2])

            if length(tokens) >= 3
                z = tryparse(T, tokens[3])

                if length(tokens) >= 4
                    warn && @warn "More than three coordinates specified in SEG2 trace header '$key' (trace number $trace_number)"
                end
            else
                warn && @warn "Two coordinates specified in SEG2 trace header '$key', not one or three (trace number $trace_number)"
                z = nothing
            end
        else
            y = z = nothing
        end
    else
        x = y = z = nothing
    end

    return _to_meters.(something.((x, y, z), missing), unit)
end
_parse_warn_location(::Type{T}, meta::Seis.SeisDict, key::Symbol, unit::AbstractString; kwargs...) where T =
    _parse_warn_location(T, meta[key], unit; key=key, kwargs...)
_parse_warn_location(::Type{T}, ::Missing, unit::AbstractString; kwargs...) where T =
    missing, missing, missing

"""
    _to_meters(value, unit; warn=true) -> value_in_m

Convert `value` into a quantity in meters according to the contents
of the string `unit`.
"""
function _to_meters(value, unit; warn=true)
    if unit == "METERS" || unit == "NONE"
        value
    elseif unit == "CENTIMETERS"
        100*value
    elseif unit == "INCH"
        254//10000*value
    elseif unit == "FEET"
        3048//10000*value
    else
        warn && @warn "Cannot interpret unit \"$unit\"; assuming distances are in m"
        value
    end
end

"Date formats to try when parsing `ACQUISITION_DATE` in the file descriptor strings"
const _TRIAL_DATEFORMATS = (
    Dates.dateformat"yyyy-mm-dd",
    Dates.dateformat"d/u/y",
    Dates.dateformat"d/U/y",
    Dates.dateformat"d/m/y",
)

"""
    _tryparse_date(fd::FileDescriptorBlock, default=nothing) -> ::Union{Nothing,Dates.DateTime}

Try to parse the date and time from the file descriptor of a SEG2 file.
The following date formats are attemped (in this order) before giving
up and returning `nothing` if none match or there is no acquisition
time specified: $(_TRIAL_DATEFORMATS)

By preference, use and explicit string labelled `"ACQUISITION_DATE_UTC"`
over `"ACQUISITION_DATE"` and likewise for the time.
"""
function _tryparse_time(fd::FileDescriptorBlock, default=nothing)
    date::Union{Nothing,Dates.Date} = nothing
    # Whether we have found
    explicit_utc = false

    if haskey(fd.strings, "ACQUISITION_DATE_UTC") || haskey(fd.strings, "ACQUISITION_DATE")
        date_string, explicit_utc = if haskey(fd.strings, "ACQUISITION_DATE_UTC")
            fd.strings["ACQUISITION_DATE_UTC"], true
        else
            fd.strings["ACQUISITION_DATE"], false
        end

        isempty(date_string) && return default
        for fmt in _TRIAL_DATEFORMATS
            maybe_date = tryparse(Dates.Date, date_string, fmt)
            if !isnothing(maybe_date)
                date = maybe_date
                break
            end
        end

        isnothing(date) && return default
    else
        return default
    end

    if haskey(fd.strings, "ACQUISITION_TIME_UTC") || haskey(fd.strings, "ACQUISITION_TIME")
        time_string = if haskey(fd.strings, "ACQUISITION_TIME_UTC")
            if !explicit_utc
                return default
            else
                fd.strings["ACQUISITION_TIME_UTC"]
            end
        else
            if explicit_utc
                return default
            else
                fd.strings["ACQUISITION_TIME"]
            end
        end

        isempty(time_string) && return default
        maybe_time = tryparse(Dates.Time, strip(time_string))
        if !isnothing(maybe_time)
            return Dates.DateTime(date, maybe_time)
        end
    else
        return default
    end
end
