"""
# `SeisMakieExt`

Extension module providing plots using Makie when a Makie backend
module is loaded.
"""
module SeisMakieExt

import Dates
import Seis

@static if isdefined(Base, :get_extension)
    import Makie
else
    import ..Makie
end

"""
    plot_traces(::AbstractArray{<:Seis.AbstractTrace}; kwargs...) -> ::Makie.Figure
    plot_traces(::AbstractTrace; kwargs...) -> ::Makie.Figure

Plot a set of traces as a set of separate wiggles, each with its own axis.

## Keyword arguments
These affect the way that the traces are plotted and annotated.

- `label=Seis.channel_code`: Set the label for each trace, which is placed
  in the top right corner of its axis.  The behaviour of this keyword depends
  on the type of `label`:
  - `::Symbol`: Take values from the `meta` dictionary if each traces
  - `::AbstractArray`: Values are taken from each entry in `label`
  - Otherwise, `label` is assumed to be a function or callable object which
    returns a string 
- `show_picks=true`: If `false`, do not plot picks.
- `sort=nothing`: Sort the traces according to one of the following:
  - `:dist`: Epicentral distance
  - `:alpha`: Alphanumerically by channel code
  - `::AbstractVector`: According to the indices in a vector passed in, e.g. from
    a call to `sortperm`.
  - Default: no sorting
- `ylims=nothing`: Control y-axis limits of traces:
  - `:all`: All traces have same amplitude limits
  - Default: each trace's axes match the trace's limits

## Makie keyword arguments
These affect the way Makie draws the figure, axes and lines.

- `figure=(size=(700,800),)`: Keyword arguments passed to the `Makie.Figure`
  constructor.
- `axis=(xgridvisible=false, ygridvisible=false)`: Keyword arguments passed
  to the `Makie.Axis` constructor.
- `lines=()`: Keyword arguments passed to the `Makie.lines` function which
  displays trace lines.
"""
function Seis.plot_traces(
    ts::AbstractArray{<:Seis.AbstractTrace};
    figure=(size=(700,800),),
    axis=(xgridvisible=false, ygridvisible=false),
    lines=(),
    label=Seis.channel_code,
    show_picks=true,
    sort=nothing,
    ylims=nothing,
    # TODO: Remove deprecated keyword arguments in a new release
    fig_kwargs=nothing,
    ax_kwargs=nothing,
    lines_kwargs=nothing,
)
    isempty(ts) && throw(ArgumentError("cannot plot empty array of traces"))

    figure = _kwargs_deprecation(:fig_kwargs, :figure, fig_kwargs, figure)
    axis = _kwargs_deprecation(:ax_kwargs, :axis, ax_kwargs, axis)
    lines = _kwargs_deprecation(:lines_kwargs, :lines, lines_kwargs, lines)

    ntraces = length(ts)

    # Sort traces into desired order for plotting
    ts = if !isnothing(sort)
        if sort == :alpha
            Base.sort(ts, by=Seis.channel_code)
        elseif sort == :dist
            Base.sort(ts, by=Seis.distance_deg)
        elseif sort isa AbstractArray
            if length(sort) != ntraces
                throw(ArgumentError(
                    "Length of `sort` ($(length(sort))) is not number of traces ($ntraces)"
                ))
            end

            ts[vec(sort)]
        else
            throw(ArgumentError("value of `sort` ($sort) not recognised"))
        end
    else
        ts
    end

    # Default axis limits
    limits = (minimum(Seis.starttime, ts), maximum(Seis.endtime, ts), nothing, nothing)
    # Custom limits
    if ylims == :all
        minval = minimum(minimum∘Seis.trace, ts)
        maxval = maximum(maximum∘Seis.trace, ts)
        limits = (limits[1], limits[2], minval, maxval)
    elseif !isnothing(ylims)
        throw(ArgumentError("valid values for `ylims` are: `:all`"))
    end

    # Get trace labels
    labels = if label isa AbstractArray
        length(label) == ntraces ||
            throw(ArgumentError("number of labels must equal the number of traces"))
        label
    elseif label isa Symbol
        [coalesce(t.meta[label], "") for t in ts]
    else
        label.(ts)
    end

    fig = Makie.Figure(; figure...)

    axs = Vector{Makie.Axis}(undef, length(ts))

    # Plot traces
    for (i, t) in enumerate(ts)
        axs[i] = Makie.Axis(fig[i,1]; limits=limits, axis...)
        Makie.lines!(axs[i], Seis.times(t), Seis.trace(t);
            color=:black, linewidth=1, lines...
        )
        Makie.linkxaxes!(axs[i], axs[1])

        if i < ntraces
            Makie.hidexdecorations!(axs[i])
        end

        # Trace labels in the top right
        if !isnothing(labels)
            Makie.text!(axs[i], 1, 1;
                text=labels[i],
                space=:relative,
                align=(:right, :top),
                offset=(-2, -1)
            )
        end

        # Picks
        if show_picks
            picks = t.picks
            if !isempty(picks)
                pick_keys = collect(keys(picks))
                pick_times = [picks[key].time for key in pick_keys]
                pick_names = String[
                    coalesce(picks[key].name, String(key)) for key in pick_keys
                ]

                Makie.vlines!(axs[i], pick_times, color=:blue)
                # Work around inability to place things at a relative position
                # in one dimension and absolute in another by plotting the
                # pick label at the bottom of the axis in data space, using
                # `Axis.finallimits`.  This does mean that if you zoom the
                # axis in the y direction, the pick label moves
                # FIXME: Find a better way to plot the label at a constant
                #        y-position.
                ys = fill(axs[i].finallimits[].origin[2], length(pick_keys))
                Makie.text!(axs[i], pick_times, ys;
                    text=pick_names, color=:blue, fontsize=8, offset=(2, 0))
            end
        end
    end

    Makie.rowgap!(fig.layout, 0)

    axs[end].xlabel = if isempty(axis)
        "Time / s"
    else
        get(axis, :xlabel, "Time / s")
    end

    fig
end

Seis.plot_traces(t::Seis.AbstractTrace; kwargs...) = Seis.plot_traces([t]; kwargs...)

function Makie.plot(ts::AbstractArray{<:Seis.AbstractTrace}; kwargs...)
    Seis.plot_traces(ts; kwargs...)
end
Makie.plot(t::Seis.AbstractTrace; kwargs...) = Makie.plot([t]; kwargs...)

"""
    plot_section!([ax::Makie.Axis=Makie.current_axis(),] traces::AbstractVector{<:Seis.AbstractTrace}, y_values=Seis.distance_deg; kwargs...)

Plot a record section for the `traces` supplied, where each is plotted at
`y_values` against time on the x-axis.  If no explicit `ax` is given,
then the most recently used `Makie.Axis` is updated.

`y_values` can be one of the following:
- A function, in which case `y_values` is applied to each trace and the value
  of that function is used.
- A `Symbol`, in which case the entry for each trace's `.meta` field with that
  key is used at the value.
- An `AbstractArray` of values, where the `i`th value of `y_values[i]` is
  used for the `i`th trace `traces[i]`.
- A string:
  - `"index"`: Plot equally spaced apart by trace index, starting at 1.

## Keyword arguments
- `absscale`: Set to a value to plot traces at some absolute scale.  This is
  useful if one wants two or more sections to have the same scale.
- `align`: Set to a `String` to align on the first pick with this name.
  Set to an array of values to align on the value for each trace.
  Set to a `Symbol` to use the pick of each trace with that key.
- `color = :black`: Line color for traces; passed to `Makie.lines`.
- `decimate`: If `false`, do not perform downsampling of traces for plotting.
  Defaults to `true`.
- `lines`: Passed to the `Makie.lines` call which plots the traces.
- `flip = false`: Flip the polarity of traces if `true`, so that positive values
  point down the page.  Note that positive values are always up even if
  `reverse` is `true`, unless `flip` is also `true`.
- `fill_down`: Set the fill colour of the negative parts of traces.  Passed
  as the `color` keyword argument to `Makie.band`.
- `fill_up`: Set the fill colour of the positive parts of traces.  Passed
  as the `color` keyword argument to `Makie.band`.
  (Use `linewidth=0` to turn off drawing of lines with the `fill` options.)
- `linewidth = 1`: Width of trace lines.  Passed to `Makie.lines`.
  If `linewidth` is ≤ 0, no traces are plotted.
- `max_samples = 1_000_000`: Control the maximum number of samples to display
  at one time in order to make plotting quicker.  Set `decimate` to `false` to
  turn this off.
- `reverse = false`: If `true`, reverse the sense of the y axis such that values
  increase down the plot rather than up.
- `show_picks`:  If `true`, add marks on the record section for each pick in the trace
  headers.
- `zoom`: Set magnification scale for traces (default 1).

## Interactions
The following keys can be used in addition to the standard Makie
plot interactions (when using an interactive backend):

- `-` key: reduce the trace amplitude by a constant factor.
- `=` key: Increase the trace amplitude by a constant factor.

See also: [`plot_section`](@ref plot_section).
"""
function Seis.plot_section!(
    ax::Makie.Axis,
    ts::AbstractVector{<:Seis.AbstractTrace},
    y_values=Seis.distance_deg;
    absscale=nothing,
    align=nothing,
    color=:black,
    decimate=true,
    # fill_down=nothing,
    # fill_up=nothing,
    flip=false,
    lines=(),
    linewidth=1,
    max_samples=50_000,
    show_picks=false,
    zoom=1,
    # TODO: Remove deprecated keyword arguments
    lines_kwargs=nothing,
)
    isempty(ts) && throw(ArgumentError("set of traces cannot be empty"))

    lines = _kwargs_deprecation(:lines_kwargs, :lines, lines_kwargs, lines)

    ntraces = length(ts)
    total_npts = sum(Seis.nsamples, ts)

    shifts = Seis.Plot.time_shifts(ts, align)

    y_shifts = _calculate_y_shifts(ts, y_values)

    # Sort so bottom traces plot last
    order = sortperm(y_shifts, rev=!ax.yreversed[])
    shifts = shifts[order]
    y_shifts = y_shifts[order]
    ts = view(ts, order)

    # Scale
    polarity = (flip ⊻ ax.yreversed[]) ? -1 : 1
    scale = isnothing(absscale) ? abs(maximum(y_shifts) - minimum(y_shifts))/10 : absscale
    scale = zoom*polarity*scale
    
    # Traces
    maxval = isnothing(absscale) ? maximum(t -> maximum(abs, Seis.trace(t)), ts) : 1
    # all_times = [Seis.times(tt) .- shift for (tt, shift) in zip(ts, shifts)]
    # time = [Seis.times(tt)[1:ndecimate:end] .- shift for (tt,shift) in zip(ts, shifts)]
    # traces = [(scale.*Seis.trace(tt)./maxval .+ y)[1:ndecimate:end] for (tt, y) in zip(ts, y_shifts)]

    # Filled portions
    # if !isnothing(fill_up) || !isnothing(fill_down)
    #     for (tt, yy, level) in zip(time, traces, y_shifts)
    #         t⁻, y⁻, t⁺, y⁺ = Seis.Plot._below_above(tt, yy, level)

    #         if !isnothing(fill_down) && !isempty(t⁻)
    #             Makie.band!(ax, t⁻, y⁻, fill(level, length(t⁻)); color=fill_down)
    #         end

    #         if !isnothing(fill_up) && !isempty(t⁺)
    #             Makie.band!(ax, t⁺, y⁺, fill(level, length(t⁺)); color=fill_up)
    #         end
    #     end
    # end

    ## Interactions
    minimum_not_nan(vals) = minimum(x -> isnan(x) ? typemax(x) : x, vals)
    maximum_not_nan(vals) = maximum(x -> isnan(x) ? typemin(x) : x, vals)

    interactive_zoom_level = Makie.Observable(0.0)
    scale_factor = Makie.@lift begin
        (2.0^$interactive_zoom_level*scale/maxval)
    end

    window_tlimits = Makie.@lift begin
        limits = $(ax.finallimits)
        (limits.origin[1], limits.origin[1] + limits.widths[1])
    end

    # If all data can be plotted at once, no need to resample
    # when axis limits change
    line_data_and_text = Makie.@lift begin
        data = Makie.Point2{Float64}[]
        fill_up_data = Makie.Point2{Float64}[]
        fill_down_data = Makie.Point2{Float64}[]
        fill_up_level = Makie.Point2{Float64}[]
        fill_down_level = Makie.Point2{Float64}[]

        window_t1 = $(window_tlimits)[1]
        window_t2 = $(window_tlimits)[2]

        npts_in_window = sum(zip(ts, shifts)) do (t, shift)
            Seis.nsamples(t, window_t1 + shift, window_t2 + shift)
        end

        if npts_in_window > max_samples
            any_downsampled = false
            if linewidth > 0
                for (t, shift, y) in zip(ts, shifts, y_shifts)
                    n = _decimation_value(t, shift, window_t1, window_t2, round(Int, max_samples/ntraces))
                    n > 2 && (any_downsampled = true)
                    times_binned, traces_binned = Seis.Plot._bin_min_max(
                        t, window_t1 + shift, window_t2 + shift, n
                    )
                    append!(data, Makie.Point2{Float64}.(
                        times_binned .- shift,
                        traces_binned.*$scale_factor .+ y
                    ))
                    push!(data, Makie.Point2(NaN, NaN))
                end
                plot_text = any_downsampled ? "Downsampled " : ""
            else
                plot_text = ""
            end

        else
            for (t, shift, y) in zip(ts, shifts, y_shifts)
                i1, i2 = Seis._cut_time_indices(
                    t, window_t1 + shift, window_t2 + shift; warn=false, allowempty=true
                )

                times = Seis.times(t)[i1:i2] .- shift
                trace = @view(Seis.trace(t)[i1:i2]).*$scale_factor .+ y

                # Lines
                if linewidth > 0
                    append!(data, Makie.Point2{Float64}.(times, trace))
                    push!(data, Makie.Point2(NaN, NaN))
                end

                # Filled parts
                #=
                if !isnothing(fill_up) || !isnothing(fill_down)
                    t⁻, y⁻, t⁺, y⁺ = Seis.Plot._below_above(times, trace, y)

                    if !isnothing(fill_down) && !isempty(t⁻)
                        append!(fill_down_data, Makie.Point2{Float64}.(t⁻, y⁻))
                        push!(fill_down_data, Makie.Point2{Float64}(NaN, NaN))

                        append!(fill_down_level, Makie.Point2{Float64}.(t⁻, fill(y, length(t⁻))))
                        push!(fill_down_level, Makie.Point2{Float64}(NaN, NaN))
                    end

                    if !isnothing(fill_up) && !isempty(t⁺)
                        append!(fill_up_data, Makie.Point2{Float64}.(t⁺, y⁺))
                        push!(fill_up_data, Makie.Point2{Float64}(NaN, NaN))

                        append!(fill_up_level, Makie.Point2{Float64}.(t⁺, fill(y, length(t⁺))))
                        push!(fill_up_level, Makie.Point2{Float64}(NaN, NaN))
                    end
                end
                =#
            end


            plot_text = ""
        end

        data, fill_down_data, fill_down_level, fill_up_data, fill_up_level, plot_text
    end

    # Fills
    #=
    if !isnothing(fill_down)
        fill_down_data = Makie.@lift($(line_data_and_text)[2])
        fill_down_level = Makie.@lift($(line_data_and_text)[3])
        Makie.band!(ax, fill_down_data, fill_down_level; color=fill_down)
    end

    if !isnothing(fill_up)
        fill_up_data = Makie.@lift($(line_data_and_text)[4])
        fill_up_level = Makie.@lift($(line_data_and_text)[5])
        Makie.band!(ax, fill_up_data, fill_up_level; color=fill_up)
    end
    =#

    # Lines
    if linewidth > 0
        line_data = Makie.@lift($(line_data_and_text)[1])
        pl = Makie.lines!(ax, line_data; linewidth, color, lines...)
    end


    # Show if data are downsampled
    coarse_plot_text = Makie.@lift($(line_data_and_text)[6])
    Makie.text!(1, 0;
        text=coarse_plot_text, space=:relative, align=(:right, :bottom),
        fontsize=10
    )

    # Picks
    if show_picks
        picks = reduce(
            vcat,
            ((time - t_shift, y, coalesce(name, string(key)))
                for (tt, y, t_shift) in zip(ts, y_shifts, shifts)
                for (key, (time, name)) in tt.picks)
        )
        picks_t = first.(picks)
        picks_y = getindex.(picks, 2)
        picks_text = last.(picks)

        Makie.scatter!(ax, picks_t, picks_y; marker=:vline, markersize=14, color=:red)
        Makie.text!(
            ax, picks_t, picks_y; text=picks_text, fontsize=10,
            align=(:center, :bottom), offset=(0, 3)
        )
    end

    # Zoom traces: `=` for in, `-` for out
    Makie.on(Makie.events(ax).keyboardbutton) do event
        event.key in (Makie.Keyboard.equal, Makie.Keyboard.minus) || return
        if event.action == Makie.Keyboard.release
            is_zoom_in = event.key == Makie.Keyboard.equal
            interactive_zoom_level[] += is_zoom_in ? 1 : -1
        end
    end
end

function Seis.plot_section!(
    ts::AbstractArray{<:Seis.AbstractTrace},
    y_values=Seis.distance_deg.(ts);
    kwargs...
)
    Seis.plot_section!(Makie.current_axis(), ts, y_values; kwargs...)
end

"""
    plot_section(traces::AbstractArray{<:Seis.AbstractTrace}, y_values=Seis.distance_deg; kwargs...) -> ::Makie.Figure

Create a new figure and fill it with a record section of the `traces`.
Most `kwargs` are passed onto [`plot_section!`](@ref); see
[`plot_section!`](@ref) for more information.

The following keyword arguments are unique to this function and are
not passed on:

- `figure`: Dictionary, named tuple or set of pairs containing
  keyword arguments which are passed to `Makie.Figure`, controlling
  the appearance of the figure.
- `axis`: Dictionary, named tuple or set of pairs containing
  keyword arguments which are passed to `Makie.Axis`, controlling
  the appearance of the axis.

## Interactions
The following keys can be used in addition to the standard Makie
plot interactions (when using an interactive backend):

- `-` key: reduce the trace amplitude by a constant factor.
- `=` key: Increase the trace amplitude by a constant factor.

See also: [`plot_section!`](@ref).
"""
function Seis.plot_section(
    ts::AbstractArray{<:Seis.AbstractTrace},
    y_values=Seis.distance_deg;
    align=nothing,
    lines=(),
    reverse=false,
    figure=(size=(800, 1100),),
    axis=(xlabel="Time / s", yreversed=reverse),
    # TODO: Remove deprecated keyword arguments
    fig_kwargs=nothing,
    ax_kwargs=nothing,
    lines_kwargs=nothing,
    # Other keyword arguments passed to `plot_section!`
    kwargs...
)
    figure = _kwargs_deprecation(:fig_kwargs, :figure, fig_kwargs, figure)
    axis = _kwargs_deprecation(:ax_kwargs, :axis, ax_kwargs, axis)
    lines = _kwargs_deprecation(:lines_kwargs, :lines, lines_kwargs, lines)

    fig = Makie.Figure(; figure...)

    ylabel = if y_values isa Symbol
        String(y_values)
    elseif y_values isa AbstractString && y_values == "index"
        "Trace number"
    elseif y_values == Seis.distance_deg
        "Distance / °"
    elseif y_values == Seis.distance_km
        "Distance / km"
    else
        ""
    end

    y_shifts = _calculate_y_shifts(ts, y_values)

    # It's ugly that we have to do this both here and in the
    # mutating method, but neater than passing these values in
    # some custom way to some internal mutating function...
    shifts = Seis.Plot.time_shifts(ts, align)

    min_time = minimum(((t, shift),) -> Seis.starttime(t) - shift, zip(ts, shifts))
    max_time = maximum(((t, shift),) -> Seis.endtime(t) - shift, zip(ts, shifts))

    y_min, y_max = extrema(y_shifts)
    Δy = y_max - y_min

    limits = (min_time, max_time, y_min - Δy/20, y_max + Δy/20)

    ax = Makie.Axis(fig[1,1]; ylabel, limits, axis...)

    Seis.plot_section!(ax, ts, y_shifts; align, lines, kwargs...)

    fig
end

function _calculate_y_shifts(ts, y_values)
    ntraces = length(ts)

    if y_values isa Function
        y_values.(ts)
    elseif y_values isa Symbol
        getproperty.(ts.meta, y_values)
    elseif y_values isa AbstractArray
        length(y_values) == ntraces ||
            throw(ArgumentError(
                "length of y values ($(length(y_values))) does not " *
                "match number of traces ($ntraces)"
            ))
        y_values
    elseif y_values isa AbstractString
        # Custom
        if y_values == "index"
            1:ntraces
        else
            throw(ArgumentError("unrecognised y axis name '$y_values'"))
        end
    end
end

"""
    Seis.plot_hodogram!([ax::Makie.Axis,] t1, t2; kwargs...) -> plot

Plot the particle motion for two traces `t1` and `t2` into an existing
axis `ax` (or the currently active axis if none is given).

## Keyword arguments
- `backazimuth=false`: If `true`, plot a line showing the direction of the
  backazimuth from the first trace to the event if present in the trace
  headers.
- `backazimuth_lines=(color=:red,)`: Passed to the `Makie.lines` call
  which plots the backazimuth line.
- `lines=(color=:black,)`: Passed to the `Makie.lines` call which plots the traces.

## Example
```
julia> fig = Makie.Figure()

julia> ax = Makie.Axis(fig[1,1]; aspect=Makie.DataAspect())

julia> plot_hodogram!(ax, sample_data(:local)[1:2]...)
```

See also: [`plot_hodogram`](@ref).
"""
function Seis.plot_hodogram!(
    ax::Makie.Axis, t1::Seis.AbstractTrace, t2::Seis.AbstractTrace;
    backazimuth=false,
    lines=(color=:black,),
    backazimuth_lines=(color=:red, linewidth=2),
)
    _check_hodogram_args(t1, t2)

    maxval = _max_abs_value(t1, t2)

    pl = Makie.lines!(ax, Seis.trace(t1), Seis.trace(t2); lines...)

    if backazimuth
        _plot_hodogram_backazimuth!(ax, t1, t2, maxval, backazimuth_lines)
    end

    pl
end
Seis.plot_hodogram!(t1::Seis.AbstractTrace, t2::Seis.AbstractTrace; kwargs...) = 
    Seis.plot_hodogram!(Makie.current_axis(), t1, t2; kwargs...)

"""
    Seis.plot_hodogram!([ax::Union{Makie.Axis3,Makie.LScene},] t1, t2, t3; lines=(color=:black, linewidth=2))

Plot a 3D hodogram of three traces into an existing 3D Makie axis, if given,
or the current axis if not.  In this case the current axis must be one
of the supported Makie axis types.

# Example
```
julia> fig = Makie.Figure()

julia> ax = Makie.Axis3(fig[1,1]; aspect=(1, 1, 1))

julia> plot_hodogram!(ax, sample_data(:local)[1:3]...)
```
"""
function Seis.plot_hodogram!(
    ax::Union{Makie.Axis3,Makie.LScene},
    t1::Seis.AbstractTrace,
    t2::Seis.AbstractTrace,
    t3::Seis.AbstractTrace;
    lines=(color=:black,),
)
    _check_hodogram_args(t1, t2, t3)

    Makie.lines!(ax, Seis.trace.((t1, t2, t3))...; lines...)
end
function Seis.plot_hodogram!(
    t1::Seis.AbstractTrace, t2::Seis.AbstractTrace, t3::Seis.AbstractTrace;
    kwargs...
)
    Seis.plot_hodogram!(Makie.current_axis(), t1, t2, t3; kwargs...)
end

"""
    Seis.plot_hodogram(t1, t2; kwargs...) -> fig, ax, pl

Create a new figure and fill it with a particle motion plot or
'hodogram' of the two `Seis.AbstractTrace`s `t1` and `t2`.
Most `kwargs` are passed onto [`plot_hodogram!`](@ref); see
[`plot_hodogram!`](@ref) for more information.

The following keyword arguments are unique to this function and are
not passed on:

- `figure`: Dictionary, named tuple or set of pairs containing
  keyword arguments which are passed to `Makie.Figure`, controlling
  the appearance of the figure.
- `axis`: Dictionary, named tuple or set of pairs containing
  keyword arguments which are passed to `Makie.Axis`, controlling
  the appearance of the axis.

# Example
```
julia> ts = cut!.(sample_data(:regional)[1:2], 0, 10);

julia> plot_hodogram(ts[1], ts[2])
```

See also: [`plot_hodogram!`](@ref).
"""
function Seis.plot_hodogram(
    t1::Seis.AbstractTrace, t2::Seis.AbstractTrace;
    backazimuth=false,
    figure=(size=(320, 300),),
    axis=(aspect=Makie.DataAspect(), xgridvisible=false, ygridvisible=false),
    lines=(color=:black,),
    backazimuth_lines=(color=:red,),
)
    _check_hodogram_args(t1, t2)

    maxval = _max_abs_value(t1, t2)
    limits = 1.01*maxval.*(-1, 1, -1, 1)
    xlabel = coalesce(t1.sta.cha, string(t1.sta.azi))
    ylabel = coalesce(t2.sta.cha, string(t2.sta.azi))

    fig = Makie.Figure(; figure...)
    ax = Makie.Axis(fig[1,1];
        xlabel,
        ylabel,
        limits,
        axis...
    )

    pl = Makie.lines!(ax, Seis.trace(t1), Seis.trace(t2); lines...)

    if backazimuth
        _plot_hodogram_backazimuth!(ax, t1, t2, maxval, backazimuth_lines)
    end

    Makie.FigureAxisPlot(fig, ax, pl)
end

"""
    Seis.plot_hodogram(t1, t2, t3; figure=(size=(500, 500),), axis_type=Makie.Axis3, axis=(), lines=(color=:black, linewidth=2)) -> fig, ax, pl

Plot a 3D hodogram of three components.

As well as the `figure` and `axis` keyword arguments which are passed on
as above, this 3D version also supports the follow keyword arguments:

- `axis_type`: Makie type of 3D axis.  Currently supported are:
  - `Makie.Axis3`: A customisable 3D axis suitable for publication-quality
    plots.
  - `Makie.LScene`: A basic 3D axis which supports zooming, unlike `Axis3`.

# Example
```
julia> ts = sample_data(:regional);

julia> plot_hodogram(cut.(ts[1:3], 0, 30)...)
"""
function Seis.plot_hodogram(
    t1::Seis.AbstractTrace,
    t2::Seis.AbstractTrace,
    t3::Seis.AbstractTrace;
    figure=(size=(500, 500),),
    axis_type=Makie.Axis3,
    axis=(),
    lines=(color=:black,),
)
    _check_hodogram_args(t1, t2, t3)

    maxval = _max_abs_value(t1, t2, t3)
    limits = 1.01*maxval.*(-1, 1, -1, 1, -1, 1)
    xlabel = coalesce(t1.sta.cha, string(t1.sta.azi))
    ylabel = coalesce(t2.sta.cha, string(t2.sta.azi))
    zlabel = coalesce(t3.sta.cha, string(t3.sta.azi))

    fig = Makie.Figure(; figure...)

    axis_defaults = if axis_type == Makie.Axis3
        (; xlabel, ylabel, zlabel, limits, aspect=(1, 1, 1), viewmode=:fit)
    elseif axis_type == Makie.LScene
        ()
    else
        throw(ArgumentError("unsupported axis type $axis_type for three-component hodogram"))
    end

    ax = axis_type(fig[1,1]; axis_defaults..., axis...)

    pl = Makie.lines!(ax, Seis.trace(t1), Seis.trace(t2), Seis.trace(t3); lines...)

    Makie.FigureAxisPlot(fig, ax, pl)
end

"""
    _decimation_value(trace, shift, t1, t2, max_samples) -> n

Return the decimation value `n` required to reduce the number
of samples in `trace` between `t1` and `t2`, with a time shift
of `shift` s, to max_samples
"""
function _decimation_value(trace::Seis.AbstractTrace, shift, t1, t2, max_samples)
    npts = Seis.nsamples(trace, t1 + shift, t2 + shift)
    max(1, round(Int, npts/max_samples))
end
function _decimation_value(times, t1, t2, max_samples)
    npts = sum(x -> t1 <= x <= t2, times)
    max(1, round(Int, npts/max_samples))
end

function _check_hodogram_args(traces...)
    t1, ts = Iterators.peel(traces)

    if !all(t -> Seis.nsamples(t) == Seis.nsamples(t1), ts)
        throw(ArgumentError("both traces must be the same length"))
    elseif !all(t -> t.delta == t1.delta, ts)
        throw(ArgumentError("both traces must have the same sampling interval"))
    elseif !all(t -> Seis.starttime(t) == Seis.starttime(t1), ts)
        throw(ArgumentError("both traces must have the same start time"))
    end

    nothing
end

"""
Add the backazimuth to an existing plot, where `maxval` is the maximum absolute
amplitude across both traces.
"""
function _plot_hodogram_backazimuth!(ax, t1, t2, maxval, backazimuth_lines)
    β = Seis.backazimuth(t1) - t2.sta.azi
    xβ, yβ = (√2*maxval) .* sincos(deg2rad(β))
    Makie.lines!(ax, [0, xβ], [0, yβ]; backazimuth_lines...)
end

"Return the maximum absolute value over all traces"
_max_abs_value(ts...) = maximum(t -> maximum(abs, Seis.trace(t)), ts)

function _kwargs_deprecation(old_name, new_name, kwargs_old, kwargs_new)
    if !isnothing(kwargs_old)
        @warn("$old_name is deprecated in favour of $new_name and will be removed in a future version=")
        kwargs_old
    else
        kwargs_new
    end
end

end # module
