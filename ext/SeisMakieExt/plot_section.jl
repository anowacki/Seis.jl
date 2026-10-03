# Docstrings are `in src/makie.jl`

function Seis.plot_section!(
    ax::Makie.Axis,
    ts::AbstractVector{<:Seis.AbstractTrace},
    y_values=Seis.distance_deg;
    absscale=nothing,
    align=nothing,
    clip=false,
    color=:black,
    decimate=true,
    # TODO: Implement trace filling
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

    # Clipping
    y_clip_mins, y_clip_maxs = if clip
        min_trace_separation_y = minimum(abs, diff(sort(y_shifts)))
        y_shifts .- min_trace_separation_y/2.2, y_shifts .+ min_trace_separation_y/2.2
    else
        fill(-Inf, ntraces), fill(Inf, ntraces)
    end

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
                for (t, shift, y, y_min, y_max) in zip(ts, shifts, y_shifts, y_clip_mins, y_clip_maxs)
                    n = _decimation_value(t, shift, window_t1, window_t2, round(Int, max_samples/ntraces))
                    n > 2 && (any_downsampled = true)
                    times_binned, traces_binned = Seis.Plot._bin_min_max(
                        t, window_t1 + shift, window_t2 + shift, n
                    )
                    append!(data, Makie.Point2{Float64}.(
                        times_binned .- shift,
                        clamp.(traces_binned.*$scale_factor .+ y, y_min, y_max)
                    ))
                    push!(data, Makie.Point2(NaN, NaN))
                end
                plot_text = any_downsampled ? "Downsampled " : ""
            else
                plot_text = ""
            end

        else
            for (t, shift, y, ymin, ymax) in zip(ts, shifts, y_shifts, y_clip_mins, y_clip_maxs)
                i1, i2 = Seis._cut_time_indices(
                    t, window_t1 + shift, window_t2 + shift; warn=false, allowempty=true
                )

                times = Seis.times(t)[i1:i2] .- shift
                trace = clamp.(@view(Seis.trace(t)[i1:i2]).*$scale_factor .+ y, ymin, ymax)

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
