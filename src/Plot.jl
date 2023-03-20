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
import DSP

"Default maximum number of samples to plot at one time"
const MAX_SAMPLES = 30_000
"Default decimation behaviour.  May be useful to change if using a Plots backend which
does its own decimation or doesn't suffer performance issues."
const DECIMATE = Ref(true)

"""
    plot(traces::AbstractArray{<:Seis.Trace}; kwargs...) -> ::Plots.Plot
    plot(trace::AbstractTrace; kwargs...) -> ::Plots.Plot

Plot a set of traces as a set of separate wiggles, each with its own axes.

Additional options provided:

- `decimate`: If `false`, do not decimate traces when plotting for speed.
- `fill_down`: Set the fill colour of the negative parts of traces.
- `fill_up`: Set the fill colour of the positive parts of traces.
  (Use `line=nothing` to turn off drawing of lines with the `fill` options.)
- `label`: Set to label for each trace.  This is placed in the top right
  corner.
- `show_picks`: If `false`, do not plot picks.
- `sort`: Sort the traces according to one of the following:
  - `:dist`: Epicentral distance
  - `:alpha`: Alphanumerically by channel code
  - `::AbstractVector`: According to the indices in a vector passed in, e.g. from
    a call to `sortperm`.
- `ylims`: Control y-axis limits of traces:
  - `:all`: All traces have same amplitude limits
"""
function plot end

# Recipe defining the above
@recipe function f(t::Union{Seis.AbstractTrace,AbstractArray{<:Trace}};
        decimate=DECIMATE[],
        fill_down=nothing,
        fill_up=nothing,
        max_samples=MAX_SAMPLES,
        show_picks=true,
        sort=nothing,
    )

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
    xguide --> hcat(fill("", ntraces-1)..., "Time / s")

    # Time limits
    xlims = get!(plotattributes, :xlims,
        (minimum(first.(times.(t))), maximum(last.(times.(t)))))

    # Amplitude limits for all traces
    t′ = cut.(t, xlims...; warn=false)
    if get(plotattributes, :ylims, nothing) == :all
        xmin, xmax, ymin, ymax = traces_limits(t′)
        plotattributes[:ylims] = (ymin, ymax)
    end
    all_ylims = [get(plotattributes, :ylims, extrema(trace(tt))) for tt in t′]

    # Labels
    labels = get(plotattributes, :label, channel_code.(t′))
    if labels isa AbstractString
        labels = repeat(labels, ntraces)
    elseif labels isa Symbol
        labels = [tt.meta[labels] for tt in t′]
    elseif length(labels) != ntraces
        throw(ArgumentError(
            "`label` must be a single label or have the same number of entries as traces"))
    end

    # Annotations for each subplot
    all_annotations = Dict(i=>[] for i in 1:ntraces)

    # Ticks and xlabel only on bottom x-axis
    xticks := hcat(repeat([nothing], 1, ntraces-1), :auto)

    # Decimate
    ndecimate = decimate ? decimation_value(t′, zeros(ntraces), xlims[1], xlims[2], max_samples) : 1
    traces = [trace(tt)[1:ndecimate:end] for tt in t′]
    all_times = [times(tt)[1:ndecimate:end] for tt in t′]

    # Filled portions
    if fill_up != nothing || fill_down != nothing
        for i in eachindex(traces)
            t⁻, y⁻, t⁺, y⁺ = _below_above(all_times[i], traces[i], 0)
            if fill_down != nothing
                @series begin
                    subplot := i
                    linewidth := 0
                    fillrange := 0
                    fillcolor := fill_down
                    t⁻, y⁻
                end
            end
            if fill_up != nothing
                @series begin
                    subplot := i
                    linewidth := 0
                    fillrange := 0
                    fillcolor := fill_up
                    t⁺, y⁺
                end
            end
        end
    end

    # Plot traces
    linecolor --> :black
    linewidth --> 1
    for i in eachindex(t)
        @series begin
            subplot := i
            all_times[i], traces[i]
        end
    end

    # Picks
    if show_picks
        pick_times, pick_names = _picks_in_window(t, xlims[1], xlims[end])
        # Pick lines
        seriestype := :vline
        linecolor := :blue
        linewidth --> 1
        for i in eachindex(t)
            length(pick_times[i]) == 0 && continue
            @series begin
                subplot := i
                pick_times[i]
            end
        end

        # Picks, and pick labels for later
        seriestype := :scatter
        primary := false
        annotation_params = (8, :left, :bottom, :blue)
        for i in eachindex(t)
            length(pick_times[i]) == 0 && continue
            xs = pick_times[i]
            npicks = length(xs)
            y = all_ylims[i][1]
            @series begin
                subplot := i
                append!(all_annotations[i], [(x, y, (name, annotation_params...))
                                       for (x, name) in zip(xs, pick_names[i])])
                []
            end
        end
    end

    # Labels
    label_params = (9, :black, :top, :right)
    label_x = xlims[end] - 0.007*(xlims[end] - xlims[1])
    for i in eachindex(t)
        label_y = all_ylims[i][end]
        push!(all_annotations[i], (label_x, label_y, (labels[i], label_params...)))
    end

    # Plot all labels
    # N.B.: `annotations` only seem to be plotted for the last series in the recipe
    for i in eachindex(t′)
        isempty(all_annotations[i]) && continue
        @series begin
            primary := false
            subplot := i
            annotations := all_annotations[i]
            []
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
- `fill_down`: Set the fill colour of the negative parts of traces.
- `fill_up`: Set the fill colour of the positive parts of traces.
           (Use `line=nothing` to turn off drawing of lines with the `fill` options.)
- `max_samples`: Control the maximum number of samples to display at one time
           in order to make plotting quicker.  Set `decimate` to `false` to turn
           this off.
- `show_picks`:  If `true`, add marks on the record section for each pick in the trace
           headers.
- `zoom`: Set magnification scale for traces (default 1).
"""
section

@userplot Section

@recipe function f(sec::Section; align=nothing, decimate=DECIMATE[],
        fill_down=nothing, fill_up=nothing,
        max_samples=MAX_SAMPLES,
        show_picks=false, zoom=1.0, absscale=nothing,
    )
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
    # Sort so bottom traces plot last
    order = sortperm(y_shifts, rev=true)
    shifts = shifts[order]
    y_shifts = y_shifts[order]
    t = view(t, order)
    # Scale
    polarity = get(plotattributes, :yflip, false) ? -1 : 1
    scale = isnothing(absscale) ? abs(maximum(y_shifts) - minimum(y_shifts))/10 : absscale
    scale = zoom*polarity*scale
    # Get decimation value
    t1, t2 = get(plotattributes, :xlims, (-Inf, Inf))
    ndecimate = decimate ? decimation_value(t, shifts, t1, t2, max_samples) : 1
    # Traces
    time = [times(tt)[1:ndecimate:end] .- shift for (tt,shift) in zip(t, shifts)]
    maxval = isnothing(absscale) ? maximum([maximum(abs, trace(tt)) for tt in t]) : 1
    traces = [(scale.*trace(tt)./maxval .+ y)[1:ndecimate:end] for (tt, y) in zip(t, y_shifts)]
    # Time limits of plot
    xlims = get!(plotattributes, :xlims, (minimum(first.(time)), maximum(last.(time))))

    # Filled portions
    if fill_up != nothing || fill_down != nothing
        for (tt, yy, level) in zip(time, traces, y_shifts)
            t⁻, y⁻, t⁺, y⁺ = _below_above(tt, yy, level)
            if fill_down != nothing
                @series begin
                    primary := false
                    linewidth := 0
                    fillrange := level
                    fillcolor := fill_down
                    t⁻, y⁻
                end
            end
            if fill_up != nothing
                @series begin
                    primary := false
                    linewidth := 0
                    fillrange := level
                    fillcolor := fill_up
                    t⁺, y⁺
                end
            end
        end
    end

    # Lines
    linecolor --> :black
    linewidth --> 1
    @series begin
        label := "" # Won't show in legend when blank
        time, traces
    end

    # Picks
    ptime, py, names = Float64[], Float64[], String[]
    if show_picks
        for (i, tt) in enumerate(t)
            for (key, (time, name)) in tt.picks
                xlims[1] <= time - shifts[i] <= xlims[end] || continue
                push!(ptime, time - shifts[i])
                push!(py, y_shifts[i])
                push!(names, coalesce(name, string(key)))
            end
        end
    end
    seriestype := :scatter
    markersize := 3
    markershape := :utriangle
    markercolor := :blue
    markerstrokealpha := 0.0
    primary := false
    annotation_params = (8, :hcenter, :bottom, :blue)
    # Pick symbols
    @series begin
        ptime, py
    end
    # Pick labels
    @series begin
        annotations := [(time, y, (name, annotation_params...))
                        for (time, y, name) in zip(ptime, py, names)]
        []
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
    plot_spectrogram(trace; powscale=:linear, normalize=true, normalise=normalize, kwargs...) -> ::Plots.Plot

Plot a spectrogram calculated with [`spectrogram`](@ref spectrogram(::AbstractTrace)).

`powscale` determines how power is scaled before plotting.  It may take one of
the following values:
- `:linear`: No scaling is performed
- `:log10`: Base-10 logarithm is taken
- `dB`: Decibels relative to the maximum power

In addition, if `powscale` is a subtype of `Function`, then that function is
applied elementwise on the matrix of power values before plotting.

If `normalize` or `normalise` are `true` (the default), power is scaled such
that the maximum power is 1 before any further scaling is performed.

Other keyword arguments are passed to Plots.

See also: [`spectorgram`](@ref spectrogram(::AbstractTrace))
"""
plot_spectrogram

@userplot Plot_spectrogram

@recipe function f(plot_spec::Plot_spectrogram; normalize=nothing, normalise=nothing, powscale=:linear)
    length(plot_spec.args) == 1 && plot_spec.args[1] isa DSP.Periodograms.Spectrogram ||
        throw(ArgumentError("first argument must be a spectrogram"))
    spec = only(plot_spec.args)

    normalise = if isnothing(normalize) && isnothing(normalise)
        true
    elseif !isnothing(normalise)
        normalise
    elseif !isnothing(normalize)
        normalize
    elseif normalize != normalise
        throw(ArgumentError("cannot supply both normalise and normalize"))
    end

    spec_times = DSP.time(spec)
    spec_freqs = DSP.freq(spec)
    spec_power = let power = DSP.power(spec)
        max_power = maximum(power)
        scale = normalise ? max_power : 1

        if powscale === :linear
            power./scale
        elseif powscale === :log10
            log10.(power./scale)
        elseif powscale === :dB
            10 .* log10.(power./max_power)
        elseif powscale isa Function
            powscale.(power./scale)
        else
            throw(ArgumentError("value of powscale not recognised"))
        end
    end

    xlims --> extrema(DSP.time(spec))
    ylims --> extrema(DSP.freq(spec))

    xguide --> "Time / s"
    yguide --> "Frequency / Hz"

    @series begin
        seriestype := :heatmap
        colormap --> :turbo
        spec_times, spec_freqs, spec_power
    end
end

"""
    decimation_value(t::AbstractArray{<:Trace}, shifts, t1, t2, max_samples) -> n

Return the decimation value `n` which ensures that `max_samples` samples
are contained within the time window `t1`-`t2` s relative to the values in `shifts`.
"""
function decimation_value(t::AbstractArray{<:Trace}, shifts, t1, t2, max_samples)
    n = 1
    np = sum(nsamples.(t, shifts .+ t1, shifts .+ t2))
    while np÷n > max_samples
        n += 1
    end
    n
end

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
                [first(Seis.picks(tt, align)).time for tt in t]
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
        shifts=zeros(promote_type(eltype.(t)...), length(t)))
    all_times = [times(tt) .+ shift for (tt, shift) in zip(t, shifts)]
    tmin = minimum(first.(all_times))
    tmax = maximum(last.(all_times))
    ymin = minimum(minimum.(trace.(t)))
    ymax = maximum(maximum.(trace.(t)))
    tmin, tmax, ymin, ymax
end

"""
    _below_above(x, y, level) -> x⁻, y⁻, x⁺, y⁺

Divide the data in `y`, with times at `x`, into two series, linearly extrapolated
such that points always lie at `level` wherever `y` crosses that value.
`x⁻` and `y⁻` are the respective point times below `level` and
`x⁺` and `y⁺` are those above.

Continuous series of points below (or below) `level` are separated in
`x⁻` and `y⁻` (or `x⁺` and `y⁺`) with `NaN`.

# Example
```
julia> x = 1:3;

julia> y = [0, 0, 2];

julia> Seis.Plot._below_above(x, y, 1)
[1.0, 2.0, 2.5, NaN], [0.0, 0.0, 1.0, NaN], [2.5, 3.0], [1.0, 2.0]
```
"""
function _below_above(x, y, level)
    n = length(x)
    length(x) == n || throw(DimensionMismatch("`x` and `y` must be the same length"))
    T = float(eltype(y))
    nan = T(NaN)
    x⁺ = T[]
    y⁺ = T[]
    x⁻ = T[]
    y⁻ = T[]
    sizehint!(x⁺, 3n÷2)
    sizehint!(y⁺, 3n÷2)
    sizehint!(x⁻, 3n÷2)
    sizehint!(y⁻, 3n÷2)
    last_point = :neither
    last_x = last_y = zero(T)
    first_point = true
    for (xx, yy) in zip(x, y)
        # First point
        if first_point
            if yy > level
                push!(x⁺, xx)
                push!(y⁺, yy)
                last_point = :above
            elseif yy < level
                push!(x⁻, xx)
                push!(y⁻, yy)
                last_point = :below
            else # y == level
                last_point = :neither
            end
            first_point = false
        # Above line
        elseif yy > level
            # Last point was above line: add this point
            if last_point === :above
                push!(x⁺, xx)
                push!(y⁺, yy)
            # Last point was below, so add the crossing point
            # to both and a delimiter (NaN) for plotting
            elseif last_point === :below
                x_cross = _intercept(last_x, last_y, xx, yy, level)
                y_cross = level
                push!(x⁺, x_cross)
                push!(y⁺, y_cross)
                push!(x⁺, xx)
                push!(y⁺, yy)
                push!(x⁻, x_cross)
                push!(y⁻, y_cross)
                push!(x⁻, nan)
                push!(y⁻, nan)
            # Last point was on the line, so add the last point and this point
            elseif last_point === :neither
                push!(x⁺, last_x)
                push!(y⁺, last_y)
                push!(x⁺, xx)
                push!(y⁺, yy)
            end
            last_point = :above
        # Below line
        elseif yy < level
            if last_point === :above
                x_cross = _intercept(last_x, last_y, xx, yy, level)
                y_cross = level
                push!(x⁺, x_cross)
                push!(y⁺, y_cross)
                push!(x⁺, nan)
                push!(y⁺, nan)
                push!(x⁻, x_cross)
                push!(y⁻, y_cross)
                push!(x⁻, xx)
                push!(y⁻, yy)
            elseif last_point === :below
                push!(x⁻, xx)
                push!(y⁻, yy)
            elseif last_point === :neither
                push!(x⁻, last_x)
                push!(y⁻, last_y)
                push!(x⁻, xx)
                push!(x⁻, yy)
            end
            last_point = :below
        # On line
        else # yy == level
            # Add this point and a delimiter as this is the first point on the line
            if last_point === :above
                push!(x⁺, xx)
                push!(y⁺, yy)
                push!(x⁺, nan)
                push!(y⁺, nan)
            elseif last_point === :below
                push!(x⁻, xx)
                push!(y⁻, yy)
                push!(x⁻, nan)
                push!(y⁻, nan)
            end
            # Otherwise the last point was on the line as well, and do nothing
            last_point = :neither
        end
        last_x = T(xx)
        last_y = yy
    end
    x⁻, y⁻, x⁺, y⁺
end

"""
    _intercept(x1, y1, x2, y2, level) -> x

Return the value of `x` at which the line defined by the coordinates
`(x1, y1)` and `(x2, y2)` intercepts `y=level`.

If `y1` and `y2` do not straddle `level`, then the value returned will not
lie between `x1` and `x2`, as expected.

If `x1` is equal to `x2`, the behaviour is undefined .
"""
function _intercept(x1, y1, x2, y2, level)
    dx = x2 - x1
    dy = y2 - y1
    ∇ = dy/dx
    x1 + (level - y1)/∇
end

"""
    _picks_in_window(traces, xlims, shifts=zeros(lengt(traces))) -> times, names

Return vectors of vectors of pick `times` and `names` within the time window
between `t1` and `t2` s.  The first element of `times` and `names` gives
the times and names of picks for the first trace in `traces`, and so on.
"""
function _picks_in_window(traces, t1, t2, shifts=zeros(length(traces)))
    length(traces) == length(shifts) ||
        throw(DimensionMismatch("traces and shifts not the same length"))
    pick_times = Vector{Float64}[]
    pick_names = Vector{String}[]
    for (i, (tt, shift)) in enumerate(zip(traces, shifts))
        push!(pick_times, [])
        push!(pick_names, [])
        for (key, (time, name)) in tt.picks
            if t1 <= time - shift <= t2
                push!(pick_times[end], time - shift)
                push!(pick_names[end], coalesce(name, string(key)))
            end
        end
    end
    pick_times, pick_names
end

end # module
