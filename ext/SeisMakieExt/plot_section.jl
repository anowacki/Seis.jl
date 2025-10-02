"""
    plot_section!([ax::Makie.Axis=Makie.current_axis(),] traces::AbstractVector{<:Seis.AbstractTrace}, y_values=Seis.distance_deg; kwargs...)::Makie.Plot

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
    pl = if linewidth > 0
        line_data = Makie.@lift($(line_data_and_text)[1])
        Makie.lines!(ax, line_data; linewidth, color, lines...)
    else
        Makie.lines!(ax, [NaN32], [NaN32])
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
        # If there is only one pick, then `reduce(vcat, ...)` does not return a vector:
        # https://github.com/JuliaLang/julia/issues/34380
        picks = (picks isa AbstractArray) ? picks : [picks]
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

    pl
end

function Seis.plot_section!(
    ts::AbstractArray{<:Seis.AbstractTrace},
    y_values=Seis.distance_deg;
    kwargs...
)
    Seis.plot_section!(Makie.current_axis(), ts, y_values; kwargs...)
end

"""
    plot_section(traces::AbstractArray{<:Seis.AbstractTrace}, y_values=Seis.distance_deg; kwargs...) -> (fig, ax, pl)::Makie.FigureAxisPlot

Create a new figure and fill it with a record section of the `traces`.
Most `kwargs` are passed onto [`plot_section!`](@ref); see
[`plot_section!`](@ref) for more information.

This function returns a `Makie.FigureAxisPlot` object containing the
figure handle `fig`, axis object `ax` and plot `pl`.

The following keyword arguments are unique to this function and are
not passed on to [`plot_section!`]:

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
    figure=(),
    # TODO: Remove deprecated keyword arguments
    fig_kwargs=nothing,
    # Extra kwargs passed to `plot_section(::AbstractArray{<:Seis.AbstractTrace}, ...)`
    kwargs...
)
    figure = _kwargs_deprecation(:fig_kwargs, :figure, fig_kwargs, figure)
    fig = Makie.Figure(; size=(800, 1100), figure...)
    ax, pl = Seis.plot_section(fig[1,1], ts, y_values; kwargs...)
    Makie.FigureAxisPlot(fig, ax, pl)
end

"""
    plot_section(gridposition, traces, y_values=Seis.distance_deg; kwargs...) -> (ax, pl)::Makie.AxisPlot

Create a new record section plot at `gridposition`, which is the
position within a Makie layout.  Usually, this will be created by
indexing into a `Makie.Figure` object.

Returns a `Makie.AxisPlot` object containing the axis handle `ax`
and plot object `pl`.

# Example
```
julia> import GLMakie as Makie

julia> fig = Makie.Figure();

julia> plot_section(fig[1,1], sample_data(:array))
```
"""
function Seis.plot_section(
    gp::Union{Makie.GridPosition,Makie.GridSubposition},
    ts::AbstractArray{<:Seis.AbstractTrace},
    y_values=Seis.distance_deg;
    align=nothing,
    lines=(),
    reverse=false,
    axis=(),
    # TODO: Remove deprecated keyword arguments
    ax_kwargs=nothing,
    lines_kwargs=nothing,
    # Other keyword arguments passed to `plot_section!`
    kwargs...
)
    axis = _kwargs_deprecation(:ax_kwargs, :axis, ax_kwargs, axis)
    lines = _kwargs_deprecation(:lines_kwargs, :lines, lines_kwargs, lines)

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

    ax = Makie.Axis(gp; ylabel, limits, xlabel="Time / s", yreversed=reverse, axis...)

    pl = Seis.plot_section!(ax, ts, y_shifts; align, lines, kwargs...)

    Makie.AxisPlot(ax, pl)
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
