# Stubs for plotting with Makie when it is loaded

"""
    plot_hodogram!

Function to plot the particle motion of two (in 2D) or three (in 3D)
traces.  This can be used once you have loaded a
[Makie](https://docs.makie.org/stable) backed (e.g., via
`using GLMakie`, `using CairoMakie`, etc.).

See docstrings below after loading a Makie backend for more information.
"""
function plot_hodogram! end

"""
    plot_hodogram

Function to plot the particle motion of two (in 2D) or three (in 3D)
traces.  This can be used once you have loaded a
[Makie](https://docs.makie.org/stable) backed (e.g., via
`using GLMakie`, `using CairoMakie`, etc.).

See docstrings below after loading a Makie backend for more information.
"""
function plot_hodogram end

"""
    plot_traces

Function to plot several traces which can be used once you have
loaded a [Makie](https://docs.makie.org/stable) backend (e.g.,
via `using GLMakie`, `using CairoMakie`, etc.).

See docstrings below after loading a Makie backend for more information.
"""
function plot_traces end

"""
    plot_section!

Function to plot several traces as a function of distance or some
other variable on the same axis.  This can be used once you have loaded
a [Makie](https://docs.makie.org/stable) backend (e.g., via
`using GLMakie`, `using CairoMakie`, etc.).

See docstrings below after loading a Makie backend for more information.
"""
function plot_section! end

"""
    plot_section

Function to plot several traces as a function of distance or some
other variable on the same axis.  This can be used once you have loaded
a [Makie](https://docs.makie.org/stable) backend (e.g., via
`using GLMakie`, `using CairoMakie`, etc.).

See docstrings below after loading a Makie backend for more information.
"""
function plot_section end

"""
    pick_axis(ax::Makie.Axis) -> (; time, yvalue)

Pick times on an existing `Makie.Axis` object, as returned from
any of the `plot_` plotting functions.

Returns a vector of named tuples, each with keys `time` (the pick
time in s) and `yvalue` (the value of the independent axis).
"""
function pick_axis end
