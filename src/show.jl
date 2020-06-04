# Display of types

_dict_summary(d) = (nkeys = length(d);
    string("meta: ", nkeys, nkeys == 1 ? " key" : " keys"))

## Event
function Base.show(io::IO, e::Event{T,P}) where {T,P}
    print(io, "Event: ")
    firstitem = true
    for f in propertynames(e)
        f in (:pos, :meta) && continue
        if getproperty(e, f) !== missing
            if !firstitem
                print(io, ", ")
            end
            print(io, "$f: $(getproperty(e, f))")
            firstitem = false
        end
    end
    firstitem || print(io, ", ")
    print(io, _dict_summary(e.meta))
end

function Base.show(io::IO, ::MIME"text/plain", e::Event{T,P}) where {T,P}
    hdr_string_len = maximum(length(string.(propertynames(e))))
    indent = 4
    padded_print = (name, val) ->
        print(io, "\n", lpad(string(name), hdr_string_len + indent) * ": ", val)
    print(io, "Seis.Event{$T,$P}:")
    # Non-meta fields
    for field in propertynames(e)
        field in (:meta, :pos) && continue
        padded_print(field, getproperty(e, field))
    end
    # Meta fields
    padded_print("meta", "")
    show_dict(io, e.meta, hdr_string_len, indent)
end

## Station
function Base.show(io::IO, s::Station{T,P}) where {T,P}
    print(io, "Station: ", channel_code(s))
    for f in propertynames(s)
        f in (:net, :sta, :loc, :cha, :pos, :meta) && continue
        if getproperty(s, f) !== missing
            print(io, ", $f: $(getproperty(s, f))")
        end
    end
    print(io, ", ", _dict_summary(s.meta))
end

function Base.show(io::IO, ::MIME"text/plain", s::Station{T,P}) where {T,P}
    hdr_string_len = maximum(length(string.(propertynames(s))))
    indent = 4
    padded_print = (name, val) ->
        print(io, "\n", lpad(string(name), hdr_string_len + indent) * ": ", val)
    print(io, "Seis.Station{$T,$P}:")
    # Non-meta fields
    for field in propertynames(s)
        field in (:meta, :pos) && continue
        padded_print(field, getproperty(s, field))
    end
    # Meta fields
    padded_print("meta", "")
    show_dict(io, s.meta, hdr_string_len, indent)
end

## Pick
Base.show(io::IO, p::Pick{T}) where T =
    print(io, "Seis.Pick{$T}(time=$(p.time), name=$(repr(p.name)))")

## Trace
# List printing
function Base.show(io::IO, t::Trace)
    code = channel_code(t)
    out = "Seis.Trace($code: delta=$(t.delta), b=$(t.b), nsamples=$(nsamples(t))"
    print(io, out * ")")
end

# Printing to the REPL, for instance
function Base.show(io::IO, ::MIME"text/plain", t::Trace{T,V,P}) where {T,V,P}
    hdr_string_len = maximum(length.(string.([TRACE_FIELDS...; STATION_FIELDS...;
                                              EVENT_FIELDS...]))) + 4
    indent = 4
    # Closure to print individual fields
    padded_print = (name, val) ->
        print(io, "\n", lpad(string(name), hdr_string_len + indent) * ": ", val)
    # Header for whole output
    print(io, "Seis.Trace{$T,$V,$P}:")
    for field in (:b, :delta)
        padded_print(field, getproperty(t, field))
    end
    for (item, fields) in zip((:sta, :evt), (propertynames(t.sta), propertynames(t.evt)))
        print(io, "\n", " "^(indent÷2-1), string(typeof(getproperty(t, item))), ":")
        # Prefix component, station and event fields with their fieldnames
        # within the Trace struct, giving e.g. 'evt.lon' rather than just 'lon'
        field_prefix = string(item)
        for field in fields
            field === :pos && continue
            value = getproperty(getproperty(t, item), field)
            if !ismissing(value)
                field_name = join((field_prefix, string(field)), ".")
                if value isa Dict && length(value) > 0
                    padded_print(field_name, "")
                    kindex = 0
                    for (k, v) in value
                        kindex += 1
                        if kindex > 1
                            print(io, "\n", " "^(hdr_string_len+indent+2), k, " => ", v)
                        else
                            print(io, k, " => ", v)
                        end
                    end
                else
                    padded_print(field_name, value)
                end
            end
        end
    end
    # Extra info for trace
    print(io, "\n Trace:")
    padded_print("picks", "")
    print(length(t.picks))
    padded_print("meta", "")
    show_dict(io, t.meta, hdr_string_len, indent)
end

function show_dict(io, dict, hdr_string_len, indent)
    kindex = 0
    for (k, v) in dict
        kindex += 1
        if kindex > 1
            print(io, "\n", " "^(hdr_string_len+indent+2), k, " => ")
        else
            print(io, k, " => ")
        end
        show(io, v)
    end
end
