"""
# `SeisMakieExt`

Extension module providing plots using Makie when a Makie backend
module is loaded.
"""
module SeisMakieExt

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

- `fig_kwargs=(size=(700,800),)`: Keyword arguments passed to the `Makie.Figure`
  constructor.
- `ax_kwargs=(xgridvisible=false, ygridvisible=false)`: Keyword arguments passed
  to the `Makie.Axis` constructor.
- `lines_kwargs=()`: Keyword arguments passed to the `Makie.lines` function which
  displays trace lines.
"""
function Seis.plot_traces(
    ts::AbstractArray{<:Seis.AbstractTrace};
    fig_kwargs=(size=(700,800),),
    ax_kwargs=(xgridvisible=false, ygridvisible=false),
    lines_kwargs=(),
    label=Seis.channel_code,
    show_picks=true,
    sort=nothing,
    ylims=nothing,
)
    isempty(ts) && throw(ArgumentError("cannot plot empty array of traces"))

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

    fig = Makie.Figure(; fig_kwargs...)

    axs = Vector{Makie.Axis}(undef, length(ts))

    # Plot traces
    for (i, t) in enumerate(ts)
        axs[i] = Makie.Axis(fig[i,1]; limits=limits, ax_kwargs...)
        Makie.lines!(axs[i], Seis.times(t), Seis.trace(t);
            color=:black, linewidth=1, lines_kwargs...
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

    axs[end].xlabel = "Time / s"

    fig
end

Seis.plot_traces(t::Seis.AbstractTrace; kwargs...) = Seis.plot_traces([t]; kwargs...)

end # module
