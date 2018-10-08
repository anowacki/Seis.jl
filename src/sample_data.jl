# Routines for accessing sample data

"Directory containing sample data"
const SAMPLE_DATA_DIR = joinpath(dirname(@__FILE__()), "..", "data")
"Types of sample data available"
const SAMPLE_DATA_KINDS = (:local, :regional, :teleseism, :teleseisl, :array)

"""
    sample_data() -> ::Trace
    sample_data(kind::Symbol) -> ::Array{Trace}

Return some sample data.

With no arguments, `sample` gives one trace from a local earthquake recorded in
California.

In the second form, a set of traces is returned according to the table below:

|`kind`|Description|
|:-----|:----------|
|`:local`|Livermore Valley, CA.  9 3-component stations|
|`:regional`|Nevada.  4 3-component stations|
|`:teleseism`|**Mid-period** recording of Eureka, CA event.  4 3-c stations|
|`:teleseisl`|**Long-period** recording of Eureka, CA event.  4 3-c stations|
|`:array`|Deep Fiji event.  60 vertical stations in the UK|
"""
sample_data() = read_sac(joinpath(SAMPLE_DATA_DIR, "seis.sac"))

function sample_data(kind::Symbol)
    kind in SAMPLE_DATA_KINDS ||
        throw(ArgumentError("No sample data matching '$kind'.  " *
                            "Choose from: $(SAMPLE_DATA_KINDS)"))
    file_pattern = "*"
    dir = joinpath(SAMPLE_DATA_DIR, string(kind))
    if kind in [:teleseism, :teleseisl]
        file_pattern = "*" * dir[end:end] * ".*"
        dir = dir[1:end-1]
    end
    t = read_sac(file_pattern, dir, echo=false)
    for (tt, file) in zip(t, t.meta.file)
        tt.sta.cha = split(file, ".")[end]
    end
    t
end
