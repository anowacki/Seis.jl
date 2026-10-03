# Stubs for plotting with Makie when it is loaded

const _MAKIE_LOADING_TEXT = """
!!! note
    This function can be used once you have loaded a
    [Makie](https://docs.makie.org/stable) backed (e.g., via
    `using GLMakie`, `using CairoMakie`, etc.).
"""

"""
    plot_hodogram!([ax::Makie.Axis,] t1, t2; kwargs...) -> plot

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

$(_MAKIE_LOADING_TEXT)

See also: [`plot_hodogram`](@ref).

---

    plot_hodogram!([ax::Union{Makie.Axis3,Makie.LScene},] t1, t2, t3; lines=(color=:black, linewidth=2))

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
function plot_hodogram! end

"""
    plot_hodogram(t1, t2; kwargs...) -> (fig, ax, pl)::Makie.FigureAxisPlot

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

$(_MAKIE_LOADING_TEXT)

See also: [`plot_hodogram!`](@ref).

---

    plot_hodogram(gridposition, t1, t2; kwargs...) -> (ax, pl)::Makie.AxisPlot

Create a new hodogram plot at `gridposition`, which is the
position within a Makie layout.  Usually, this will be created by
indexing into a `Makie.Figure` object.

Returns a `Makie.AxisPlot` object containing the axis handle `ax`
and plot object `pl`.

---

    plot_hodogram(t1, t2, t3; figure=(size=(500, 500),), axis_type=Makie.Axis3, axis=(), lines=(color=:black, linewidth=2)) -> fig, ax, pl

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
function plot_hodogram end

"""
    plot_traces(::AbstractArray{<:Seis.AbstractTrace}; kwargs...) -> ::Makie.Figure
    plot_traces(::AbstractTrace; kwargs...) -> ::Makie.Figure
    Makie.plot(::Union{AbstractTrace, AbstractArray{<:AbstractTrace}; kwargs...) -> ::Makie.Figure

Plot a set of traces as a set of separate wiggles, each with its own axis,
returning the figure object created.

## Keyword arguments
These affect the way that the traces are plotted and annotated.

- `label=Seis.channel_code`: Set the label for each trace, which is placed
  in the top right corner of its axis.  The behaviour of this keyword depends
  on the type of `label`:
  - `::Symbol`: Take values from the `meta` dictionary of each trace
  - `::AbstractArray`: Values are taken from each entry in `label`
  - Otherwise, `label` is assumed to be a function or callable object which
    returns a string and takes a single trace as its argument
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
- `lines=(color=:black, linewidth=1)`: Keyword arguments passed to the
  `Makie.lines` function which displays trace lines.

$(_MAKIE_LOADING_TEXT)

---

    plot_traces(gridposition, traces; kwargs...) -> ::Vector{Makie.AxisPlot}

Plot a set of traces, each in their own axis, into a subgrid layout
of an existing `Makie.Figure`.  Usually this will be created by
indexing into an existing figure object.

# Example
Plot the east and north components of a set of data in two sets of axes
next to each other:
```
julia> fig = plot_traces(filter(is_east, sample_data(:local)); lines=(; color=:red))

julia> plot_traces(fig[1,2], filter(is_north, sample_data(:local)); lines=(; color=:blue));
```
"""
function plot_traces end

# TODO: When implemented add the following to the docstring:
# - `fill_down`: Set the fill colour of the negative parts of traces.  Passed
#  as the `color` keyword argument to `Makie.band`.
#- `fill_up`: Set the fill colour of the positive parts of traces.  Passed
#  as the `color` keyword argument to `Makie.band`.
#  (Use `linewidth=0` to turn off drawing of lines with the `fill` options.)

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
- `clip = false`: If `true`, clip traces such that they do not overlap.
- `color = :black`: Line color for traces; passed to `Makie.lines`.
- `decimate`: If `false`, do not perform downsampling of traces for plotting.
  Defaults to `true`.
- `lines`: Passed to the `Makie.lines` call which plots the traces.
- `flip = false`: Flip the polarity of traces if `true`, so that positive values
  point down the page.  Note that positive values are always up even if
  `reverse` is `true`, unless `flip` is also `true`.
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

$(_MAKIE_LOADING_TEXT)

See also: [`plot_section`](@ref plot_section).
"""
function plot_section! end

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

$(_MAKIE_LOADING_TEXT)

See also: [`plot_section!`](@ref).

---

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
function plot_section end

"""
    plot_spectrogram!(ax::Makie.Axis, spec, quantity=:power; heatmap=()) -> ::Makie.Heatmap

Plot a spectrogram `spec` as returned by [`spectrogram`] as a heatmap.
Specity `quantity` as `:power` (the default) for spectral power, or
`:amplitude`.

Extra keyword arguments contained in `heatmap` are passed to the
call to `Makie.heatmap`, thus can be used e.g. to set the colour map
used for the heatmap, colour range, etc.

$(_MAKIE_LOADING_TEXT)

See also: [`plot_spectrogram`](@ref).
"""
function plot_spectrogram! end

"""
    plot_spectrogram(trace, quantity=:power; show_colorbar=true, figure, heatmap, colorbar, kwargs...) -> fig::Makie.Figure

Calculate the spectrogram for a single trace and plot both the 
trace timeseries on top, and a heatmap showing the spectrogram
below, returning the figure object.

`kwargs` are passed to the [`spectrogram`](@ref) function to create
the spectrogram.

By default the colour bar for the heatmap is shown on the right;
pass `show_colorbar=false` to disable this.

Keyowrd arguments `figure`, `axis`, `heatmap` and `colorbar` are
respectively passed to calls to `Makie.Figure`, `Makie.Axis`,
`Makie.heatmap` and `Makie.Colorbar`.  Pass named tuples of keyword
arguments which can therefore control the appearance of the plot.

$(_MAKIE_LOADING_TEXT)

See also: [`plot_spectrogram!`](@ref)

# Example
```
julia> import DSP

julia> plot_spectrogram(t, :amplitude; heatmap=(colormap=:viridis,), overlap=0.99, pad=10, window=DSP.hanning, length=0.2)
```

---

    plot_spectrogram(gridposition, spec, quantity=:power; axis, heatmap) -> (ax, pl)::Makie.AxisPlot

Plot a spectrogram `spec` as returned by [`spectrogram`] as a heatmap.
With this method, create a new axis at `gridposition`, which is the
position within a Makie layout.  Usually, this will be created by
indexing into a `Makie.Figure` object like `fig[1,1]`.

Specity `quantity` as `:power` (the default) for spectral power, or
`:amplitude`.

Keywords `axis` and `heatmap` are passed on respectively to the
calls to `Makie.Axis`, to create the axis, and `Makie.heatmap`
to create the plot.

Returns a `Makie.AxisPlot` object containing the axis handle `ax`
and plot object `pl`.

---

    plot_spectrogram(spec, quantity=:power; show_colorbar=true, figure, axis, heatmap, colorbar) -> (fig, ax, pl)::Makie.FigureAxisPlot

Plot a spectrogram by creating a new figure and axis.  By default the
colour bar is shown giving the colour scale for the heatmap; disable this
with `show_colorbar=false`.

---

    plot_spectrogram(spec, trace, quantity=:power; show_colorbar=true, figure, heatmap, colorbar) -> fig::Makie.Figure
    plot_spectrogram(trace, spec, quantity=:power; show_colorbar=true, figure, heatmap, colorbar) -> fig::Makie.Figure

Pass both a spectrogram `spec` and the `trace` from which it was
calculated to create a heatmap of the spectrogram and also the time
series above.  The first two arguments can be in either order.

"""
function plot_spectrogram end
