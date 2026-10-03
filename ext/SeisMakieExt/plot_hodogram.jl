# Docstrings are `in src/makie.jl`

# Plot into existing axis, 2D
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

# Plot into existing axis, 3D
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

# Create new figure and axis, 2D
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

# Plot into existing figure at a grid position, creating a new axis (2D)
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

# Create a new figure and axis, 3D
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

# Plot into existing figure at a grid position, creating a new axis (3D)
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
