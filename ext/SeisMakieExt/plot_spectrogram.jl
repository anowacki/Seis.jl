# Docstrings are `in src/makie.jl`

function Seis.plot_spectrogram!(
    ax::Makie.Axis,
    s::DSP.Periodograms.Spectrogram,
    quantity::Symbol=:power;
    heatmap=(),
)
    values = if quantity === :power
        s.power'
    elseif quantity === :amplitude
        sqrt.(s.power')
    else
        throw(ArgumentError("quantity must be one of `:power` or `:amplitude`"))
    end

    Makie.heatmap!(ax, s.time, s.freq, values; heatmap...)
end

function Seis.plot_spectrogram(
    gp::Union{Makie.GridPosition,Makie.GridSubposition},
    s::DSP.Periodograms.Spectrogram,
    quantity::Symbol=:power;
    axis=(),
    heatmap=(),
)
    ax = Makie.Axis(gp; xlabel="Time / s", ylabel="Frequency / Hz", axis...)
    hm = Seis.plot_spectrogram!(ax, s, quantity; heatmap)

    Makie.AxisPlot(ax, hm)
end

function Seis.plot_spectrogram(
    s::DSP.Periodograms.Spectrogram,
    quantity::Symbol=:power;
    figure=(),
    axis=(),
    heatmap=(),
    show_colorbar=true,
    colorbar=(),
)
    fig = Makie.Figure(; figure...)
    ax, hm = Seis.plot_spectrogram(fig[1,1], s, quantity; axis, heatmap)
    if show_colorbar
        label = if quantity === :power
            "Power"
        elseif quantity === :amplitude
            "Amplitude"
        end
        cb = Makie.Colorbar(fig[1,2], hm; label, colorbar...)
    end

    Makie.FigureAxisPlot(fig, ax, hm)
end

function Seis.plot_spectrogram(
    s::DSP.Periodograms.Spectrogram,
    t::Seis.AbstractTrace,
    quantity::Symbol=:power;
    figure=(),
    heatmap=(),
    show_colorbar=true,
    colorbar=(),
)
    fig, ax_hm, hm = Seis.plot_spectrogram(s, quantity; figure, heatmap, show_colorbar, colorbar)
    ((ax_trace, pl_trace),) = Seis.plot_traces(fig[0,1], [t])
    ax_trace.xaxisposition = :top
    ax_trace.limits = (extrema(s.time), nothing)

    Makie.rowsize!(fig.layout, 1, Makie.Relative(0.8))
    Makie.linkxaxes!(ax_trace, ax_hm)

    fig
end

function Seis.plot_spectrogram(
    t::Seis.AbstractTrace, s::DSP.Periodograms.Spectrogram, quantity::Symbol=:power;
    kwargs...
)
    Seis.plot_spectrogram(s, t, quantity; kwargs...)
end

function Seis.plot_spectrogram(
    t::Seis.AbstractTrace,
    quantity::Symbol=:power;
    figure=(), heatmap=(), show_colorbar=true, colorbar=(), kwargs...
)
    s = Seis.spectrogram(t; kwargs...)
    Seis.plot_spectrogram(s, t, quantity; figure, heatmap, show_colorbar, colorbar)
end
