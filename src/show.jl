# Display of types

## Event
function Base.show(io::IO, e::Event{T,S}) where {T,S}
    print(io, "Event{$T,$S}: ")
    firstitem = true
    for f in fieldnames(Event)
        if getfield(e, f) !== missing
            if !firstitem
                print(io, ", ")
            end
            print(io, "$f: $(getfield(e, f))")
            firstitem = false
        end
    end
end

function Base.show(io::IO, ::MIME"text/plain", e::Event{T,S}) where {T,S}
    hdr_string_len = maximum(length.(EVENT_FIILDS)) + 4
    indent = 4
end

## Trace
# List printing
function Base.show(io::IO, t::Trace)
    code = channel_code(t)
    out = "Seis.Trace($code: delta=$(t.delta), b=$(t.b), nsamples=$(nsamples(t))"
    print(io, out * ")")
end

# Printing to the REPL, for instance
function Base.show(io::IO, ::MIME"text/plain", t::Trace)
    hdr_string_len = maximum(length.(string.([TRACE_FIELDS...; STATION_FIELDS...;
                                              EVENT_FIELDS...]))) + 4
    indent = 4
    # Closure to print individual fields
    padded_print = (name, val) ->
        print(io, "\n", lpad(string(name), hdr_string_len + indent) * ": ", val)
    # Header for whole output
    print(io, "Seis.Trace:")
    for field in (:b, :delta)
        padded_print(field, getfield(t, field))
    end
    for (item, fields) in zip((:sta, :evt), (STATION_FIELDS, EVENT_FIELDS))
        print(io, "\n", " "^(indent÷2-1), string(typeof(getfield(t, item))), ":")
        # Prefix component, station and event fields with their fieldnames
        # within the Trace struct, giving e.g. 'evt.lon' rather than just 'lon'
        field_prefix = string(item)
        for field in fields
            value = getfield(getfield(t, item), field)
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
    padded_print("meta", "")
    kindex = 0
    for (k, v) in t.meta
        kindex += 1
        if kindex > 1
            print(io, "\n", " "^(hdr_string_len+indent+2), k, " => ", v)
        else
            print(io, k, " => ", v)
        end
    end
end