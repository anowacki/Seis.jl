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
    Seis.plot_hodogram(t1, t2; kwargs...) -> (fig, ax, pl)::Makie.FigureAxisPlot

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
julia> ts = cut!.(sample_data(:regional)[1:2], 50, 60);

julia> fig, axis, hod = plot_hodogram(ts[1], ts[2])
```

See also: [`plot_hodogram!`](@ref).
"""
function Seis.plot_hodogram(
    t1::Seis.AbstractTrace,
    t2::Seis.AbstractTrace;
    figure=(),
    kwargs...
)
    fig = Makie.Figure(; size=(320, 300), figure...)
    ax, pl = Seis.plot_hodogram(fig[1,1], t1, t2; kwargs...)
    Makie.FigureAxisPlot(fig, ax, pl)
end

"""
    plot_hodogram(gridposition, t1, t2; kwargs...) -> (ax, pl)::Makie.AxisPlot

Create a new hodogram plot at `gridposition`, which is the
position within a Makie layout.  Usually, this will be created by
indexing into a `Makie.Figure` object.

Returns a `Makie.AxisPlot` object containing the axis handle `ax`
and plot object `pl`.
"""
function Seis.plot_hodogram(
    gp::Union{Makie.GridPosition,Makie.GridSubposition},
    t1::Seis.AbstractTrace,
    t2::Seis.AbstractTrace;
    backazimuth=false,
    axis=(),
    lines=(),
    backazimuth_lines=(),
)
    _check_hodogram_args(t1, t2)

    maxval = _max_abs_value(t1, t2)
    limits = 1.01*maxval.*(-1, 1, -1, 1)
    xlabel = coalesce(t1.sta.cha, string(t1.sta.azi), t1.sta.sta)
    ylabel = coalesce(t2.sta.cha, string(t2.sta.azi), t2.sta.sta)

    ax = Makie.Axis(gp[1,1];
        aspect=Makie.DataAspect(),
        limits,
        xgridvisible=false,
        xlabel,
        ygridvisible=false,
        ylabel,
        axis...
    )

    pl = Makie.lines!(ax, Seis.trace(t1), Seis.trace(t2); color=:black, lines...)

    if backazimuth
        _plot_hodogram_backazimuth!(ax, t1, t2, maxval, (color=:red, backazimuth_lines...))
    end

    Makie.AxisPlot(ax, pl)
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
```
"""
function Seis.plot_hodogram(
    t1::Seis.AbstractTrace,
    t2::Seis.AbstractTrace,
    t3::Seis.AbstractTrace;
    figure=(),
    kwargs...
)
    fig = Makie.Figure(; size=(500, 500), figure...)
    ax, pl = Seis.plot_hodogram(fig[1,1], t1, t2, t3; kwargs...)
    Makie.FigureAxisPlot(fig, ax, pl)
end

function Seis.plot_hodogram(
    gp::Union{Makie.GridPosition,Makie.GridSubposition},
    t1::Seis.AbstractTrace,
    t2::Seis.AbstractTrace,
    t3::Seis.AbstractTrace;
    axis_type=Makie.Axis3,
    axis=(),
    lines=(),
)
    _check_hodogram_args(t1, t2, t3)

    maxval = _max_abs_value(t1, t2, t3)
    limits = 1.01*maxval.*(-1, 1, -1, 1, -1, 1)
    xlabel = coalesce(t1.sta.cha, string(t1.sta.azi), t1.sta.sta)
    ylabel = coalesce(t2.sta.cha, string(t2.sta.azi), t2.sta.sta)
    zlabel = coalesce(t3.sta.cha, string(t3.sta.azi), t3.sta.sta)

    axis_defaults = if axis_type == Makie.Axis3
        (; xlabel, ylabel, zlabel, limits, aspect=(1, 1, 1), viewmode=:fit)
    elseif axis_type == Makie.LScene
        ()
    else
        throw(ArgumentError("unsupported axis type $axis_type for three-component hodogram"))
    end

    ax = axis_type(gp[1,1]; axis_defaults..., axis...)

    pl = Makie.lines!(ax, Seis.trace(t1), Seis.trace(t2), Seis.trace(t3); color=:black, lines...)

    Makie.AxisPlot(ax, pl)
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
