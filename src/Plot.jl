"""
# Seis.Plot

Create plots of seismic data using [`Plots`](https://github.com/JuliaPlots/Plots.jl).

To use the plots enabled here, first load `Plots` with `using Plots`.  The `plot`
command to show individual traces will be enabled whether or not you have loaded
`Seis.Plot`, but the other commands can only be done when `using Seis.Plot`.

## Plot types

### `plot`
Calling `plot` with either a single `Seis.Trace` or array of traces will produce
a set of plots, one for each trace.

    plot(traces::AbstractArray{<:Seis.Trace}; kwargs...) -> ::Plots.Plot
    plot(trace::AbstractTrace; kwargs...) -> ::Plots.Plot

Plot a set of traces as a set of separate wiggles.

Additional options provided:

- `ylims`=`:all`: All traces have same amplitude limits
- `picks`: If `false`, do not plot picks.

### `section`
This produces a record section of an array of traces.  See documentation
for `section` for full details.

### `hodogram`
Particle motion for pairs of traces.
"""
module Plot

using ..Seis
using RecipesBase

"Default maximum number of samples to plot at one time"
const MAX_SAMPLES = 30_000
"Default decimation behaviour.  May be useful to change if using a Plots backend which
does its own decimation or doesn't suffer performance issues."
const DECIMATE = Ref(true)

"""
    plot(traces::AbstractArray{<:Seis.Trace}; kwargs...) -> ::Plots.Plot
    plot(trace::AbstractTrace; kwargs...) -> ::Plots.Plot

Plot a set of traces as a set of separate wiggles.

Additional options provided:

- `ylims`=`:all`: All traces have same amplitude limits
- `label`: Set to label for each trace.
- `pick`: If `false`, do not plot picks.
- `sort`: Sort the traces according to one of the following:
  - `:dist`: Epicentral distance
  - `:alpha`: Alphanumerically by channel code
  - `::AbstractVector`: According to the indices in a vector passed in, e.g. from
    a call to `sortperm`.
- `decimate`: If `false`, do not decimate traces when plotting for speed.
"""
function plot end

# Recipe defining the above
@recipe function f(t::Union{Seis.AbstractTrace,AbstractArray{<:Trace}};
        pick=true, sort=nothing, decimate=DECIMATE[], max_samples=MAX_SAMPLES)

    # Make single trace into a vector
    t isa AbstractArray || (t = [t])
    ntraces = length(t)

    # Trace ordering
    sort_indices = if sort isa Symbol
        if sort == :alpha
            sortperm(channel_code.(t))
        elseif sort == :dist
            sortperm(distance_deg.(t))
        else
            throw(ArgumentError("Value of `sort` ($sort) not recognised"))
        end
    elseif sort isa AbstractArray
        length(sort) == ntraces || throw(ArgumentError("Length of `sort` " *
            "($(length(sort))) is not number of traces ($ntraces)"))
        vec(sort)
    elseif sort === nothing
        1:ntraces
    else
        throw(ArgumentError("Value of `sort` ($sort) not recognised"))
    end
    t = t[sort_indices] # Traces now in correct order

    # Plotting defaults
    layout --> (ntraces, 1)
    legend --> false
    framestyle --> :box
    grid --> false
    linecolor --> :black
    linewidth --> 1

    # Set amplitude limits
    if get(plotattributes, :ylims, nothing) == :all
        plotattributes[:ylims] = extrema(vcat(trace.(t)...))
    end

    # Time limits
    xlims = get!(plotattributes, :xlims,
        (minimum(first.(times.(t))), maximum(last.(times.(t)))))

    # Amplitude limits for all traces
    t′ = cut.(t, xlims...; warn=false)
    if get(plotattributes, :ylims, nothing) == :all
        xmin, xmax, ymin, ymax = trace_limits(t′)
        plotattributes[:ylims] = (ymin, ymax)
    end
    all_ylims = [get(plotattributes, :ylims, extrema(trace(tt))) for tt in t′]

    # Ticks and xlabel only on bottom x-axis
    xticks := hcat(repeat([nothing], 1, ntraces-1), :auto)

    # Decimate
    ndecimate = decimate ? decimation_value(t′, zeros(ntraces), xlims[1], xlims[2], max_samples) : 1
    traces = [trace(tt)[1:ndecimate:end] for tt in t′]
    all_times = [times(tt)[1:ndecimate:end] for tt in t′]

    # Plot traces
    for i in eachindex(t)
        @series begin
            subplot := i
            all_times[i], traces[i]
        end
    end

    # Labels
    get!(plotattributes, :label, channel_code.(t))
    annot_params = (9, :black, :top, :right)
    markerstrokealpha := 0.0
    markeralpha := 0.0
    for i in eachindex(t)
        @series begin
            seriestype := :scatter
            subplot := i
            series_annotations := [Main.Plots.text(plotattributes[:label][i], annot_params...)]
            [xlims[end] - 0.007(xlims[end] - xlims[1])], [all_ylims[i][end]]
        end
    end

    # Picks
    if pick
        # Pick lines
        seriestype := :vline
        linecolor := :blue
        linewidth --> 1
        annot_params = (10, :blue, :left, :bottom)
        for i in eachindex(t)
            length(picks(t[i])) == 0 && continue
            @series begin
                guidefontcolor := :blue
                subplot := i
                [p.time for p in Seis.picks(t[i]) if xlims[1] <= p.time <= xlims[2]]
            end
        end

        # Pick labels: workaround (see https://discourse.julialang.org/t/annotations-and-line-widths-in-plots-jl-heatmaps/4259/3)
        seriestype := :scatter
        markerstrokealpha := 0.0
        # seriesalpha := 0.0
        markeralpha := 0.0
        primary := false
        annotation_params = (8, :left, :bottom, :blue)
        for i in eachindex(t)
            length(picks(t[i])) == 0 && continue
            @series begin
                subplot := i
                # FIXME: Update to a better way of implementing this: Main.Plots may break
                series_annotations := [Main.Plots.text.(coalesce(p.name, ""), annotation_params...)
                                       for p in Seis.picks(t[i])
                                       if xlims[1] <= p.time <= xlims[2]]
                x = [p.time for p in Seis.picks(t[i]) if xlims[1] <= p.time <= xlims[2]]
                y = get(plotattributes, :ylims, extrema(trace(t[i])))[1]
                x, repeat([y], length(x))
            end
        end
    end
end

"""
    section(traces, y_values=distance_deg.(traces); kwargs...) -> ::Plots.plot

Plot a record section for a collection of `Trace`s.  By default, the y-axis is
epicentral distance in degrees, but specify one of the following for `y_values`
to change this:

- A `Function` which will be applied to the traces, e.g., `distance_km`.
- A `Symbol` which corresponds to a field of the `Trace`s' `meta` dictionary.
- A set of values (i.e., `AbstractArray`)
- One of the following strings:
  - "index" for trace index (`i:length(traces)`)

Additional options provided via keyword arguments:

- `absscale`: Set to a value to plot traces at some absolute scale.  This is
              useful if one wants two or more sections to have the same scale.
- `align`: Set to a `String` to align on the first pick with this name.
           Set to an array of values to align on the value for each trace.
           Set to a `Symbol` to use the pick of each trace with that key.
- `decimate`: If `false`, do not perform downsampling of traces for plotting.
           Defaults to `true`.
- `max_samples`: Control the maximum number of samples to display at one time
           in order to make plotting quicker.  Set `decimate` to `false` to turn
           this off.
- `pick`:  If `true`, add marks on the record section for each pick in the trace
           headers.
- `zoom`: Set magnification scale for traces (default 1).
"""
section

@userplot Section

@recipe function f(sec::Section; align=nothing, decimate=DECIMATE[], max_samples=MAX_SAMPLES,
        pick=false, zoom=1.0, absscale=nothing)
    # Arguments
    t = sec.args[1]
    y_values = if length(sec.args) >= 2
        sec.args[2]
    else
        yguide --> "Distance / °"
        distance_deg.(t)
    end
    t isa AbstractArray{<:Trace} ||
        throw(ArgumentError("section requires an array of `Seis.Trace`s"))
    decimate isa Bool || throw(ArgumentError("decimate must be `true` or `false`"))
    # Defaults
    linecolor --> :black
    linewidth --> 1
    framestyle --> :box
    grid --> false
    xguide --> "Time / s"
    # Time shifts
    shifts = time_shifts(t, align)
    # Y-axis
    y_shifts = if y_values isa Function
        y_values.(t)
    elseif y_values isa Symbol
        getproperty.(t.meta, y_values)
    elseif y_values isa AbstractArray
        length(y_values) == length(t) ||
            throw(ArgumentError("Length of y values ($(length(y_values))) does not " *
                                "match number of traces ($(length(t)))"))
        y_values
    elseif y_values isa AbstractString
        # Custom
        if y_values == "index"
            1:length(t)
        else
            throw(ArgumentError("Unrecognised y axis name '$y_values'"))
        end
    end
    # Scale
    scale = isnothing(absscale) ? abs(maximum(y_shifts) - minimum(y_shifts))/10 : absscale
    scale = zoom*scale
    # Get decimation value
    t1, t2 = get(plotattributes, :xlims, (-Inf, Inf))
    ndecimate = decimate ? decimation_value(t, shifts, t1, t2, max_samples) : 1
    # Traces
    time = [times(tt)[1:ndecimate:end] .- shift for (tt,shift) in zip(t, shifts)]
    maxval = isnothing(absscale) ? maximum([maximum(abs, trace(tt)) for tt in t]) : 1
    traces = [(scale.*trace(tt)./maxval .+ y)[1:ndecimate:end] for (tt, y) in zip(t, y_shifts)]
    # Time limits of plot
    xlims = get!(plotattributes, :xlims, (minimum(first.(time)), maximum(last.(time))))

    # Plot
    @series begin
        label := "" # Won't show in legend when blank
        time, traces
    end

    # Picks
    ptime, py, names = Float64[], Float64[], String[]
    if pick
        for (i, tt) in enumerate(t)
            ps = picks(tt)
            for (time, name) in ps
                xlims[1] <= time - shifts[i] <= xlims[end] || continue
                push!(ptime, time - shifts[i])
                push!(py, y_shifts[i])
                push!(names, coalesce(name, ""))
            end
        end
    end
    seriestype := :scatter
    markersize := 3
    markershape := :utriangle
    markercolor := :blue
    markerstrokealpha := 0.0
    primary := false
    annotation_params = (8, :center, :bottom, :blue)
    @series begin
        series_annotations := Main.Plots.text.(names, annotation_params...)
        ptime, py
    end
end

"""
    hodogram(t1::Trace, t2::Trace; backazimuth=false) -> ::Plots.Plot

Plot the particle motion for a pair of traces `t1` and `t2`.

The two traces must have the same start time, length and sampling interval.

If `backazimuth` is `true`, plot the direction of the minor arc towards the event.
Requires event and station coordinates to be present in headers.
"""
hodogram

@userplot Hodogram

@recipe function f(hodogram::Hodogram; backazimuth=false)
    length(hodogram.args) == 2 || throw(ArgumentError("two traces are required as arguments"))
    all(isa.(hodogram.args, AbstractTrace)) ||
        throw(ArgumentError("arguments to hodogram must be Seis.AbstractTraces"))
    t1, t2 = hodogram.args
    aspect_ratio := :equal
    maxval = max(maximum(abs, trace(t1)), maximum(abs, trace(t2)))
    xlims --> 1.01.*(-maxval, maxval)
    ylims --> 1.01.*(-maxval, maxval)
    framestyle --> :box
    xguide --> coalesce(t1.sta.cha, string(t1.sta.azi))
    yguide --> coalesce(t2.sta.cha, string(t2.sta.azi))
    linecolor --> :black
    linewidth --> 1
    grid --> false
    ticks --> 3
    @series begin
        label --> ""
        trace(t1), trace(t2)
    end
    if backazimuth
        @series begin
            β = Seis.backazimuth(t1) - t2.sta.azi
            xβ, yβ = (√2*maxval) .* sincos(deg2rad(β))
            label --> "Backazimuth"
            linecolor --> :red
            [0, xβ], [0, yβ]
        end
    end
end

"""
    decimation_value(t::AbstractArray{<:Trace}, shifts, t1, t2, max_samples) -> n

Return the decimation value `n` which ensures that `max_samples` samples
are contained within the time window `t1`-`t2` s relative to the values in `shifts`.
"""
function decimation_value(t::AbstractArray{<:Trace}, shifts, t1, t2, max_samples)
    n = 1
    np = sum(num_points.(t, shifts .+ t1, shifts .+ t2))
    while np÷n > max_samples
        n += 1
    end
    n
end

"""
    num_points(t::Trace, t1, t2) -> np

Return the number of points `np` between `t1` and `t2` seconds relative to the
origin time of the trace `t`.
"""
num_points(t::Trace, t1, t2) = sum(t1 .<= times(t) .<= t2)

"""
    time_shifts(t, align) -> shifts

Return a vector `shifts` of time shifts in s which align the traces `t`.

When `align` is nothing, traces are not shifted.  Otherwise, the shift may be given
explicitly when `align` a vector, or when it is a `String` or `Regex`, the first
pick matching `align` will determine the shift.
"""
function time_shifts(t::AbstractArray{<:AbstractTrace}, align)
        if align == nothing
            zeros(length(t))
        elseif align isa String || align isa Regex
            shifts = try
                [first(picks(tt, align)).time for tt in t]
            catch err
                error("Error finding pick '$align' for every trace.  (Error: $err)")
            end
            shifts
        elseif align isa Symbol
            time_picks = getproperty(t.picks, align)
            any(x->x===missing, time_picks) && error("not every trace has a pick called '$align'")
            time_picks.time
        elseif align isa AbstractVector
            length(align) == length(t) ||
                throw(ArgumentError("`align` must have same length as number of traces"))
            align
        else
            throw(ArgumentError("`align` must be a pick name or vector of shifts"))
        end
end

"""
    traces_limits(t::AbstractArray{<:Trace}, shifts=zeros(length(t))) -> tmin, tmax, ymin, ymax

Return the extreme values of time and amplitude of all of the traces `t`, in an
inefficient manner.
"""
function traces_limits(t::AbstractArray{<:Trace},
        shifts=zeros(promote_type(eltype.(t.t)...), length(t)))
    all_times = [times(tt) .+ shift for (tt, shift) in zip(t, shifts)]
    tmin = minimum(first.(all_times))
    tmax = maximum(last.(all_times))
    ymin = minimum(minimum.(trace.(t)))
    ymax = maximum(maximum.(trace.(t)))
    tmin, tmax, ymin, ymax
end

end # module
