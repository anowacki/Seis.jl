# Plotting with Makie.jl

## Introduction
Seis comes with two kinds of plotting functionality.
The second, experimental kind relies on
[Makie](https://docs.makie.org/stable/).  To create plots with
Seis and Makie.jl, first install one of the Makie backends (see below):
```
import Pkg
Pkg.add("GLMakie")
```

Whenever you wish to create plots, do e.g.
```
julia> using Seis, GLMakie
```

The other Plots-based function is described [here](@ref Plotting-with-Plots.jl).

## Makie backends
The Makie plotting ecosystem relies on the user loading one or more different
'backend' packages which implement the plotting commands in different ways.  The
[Makie backend documentation](https://docs.makie.org/stable/explanations/backends/backends)
describes the different choices and how to switch between them.  For the
purposes of this documentation, however, we will consider two backends:

- GLMakie: Allows interactive plotting; can only save bitmap figures
- CairoMakie: Non-interactive; can save vector graphics for publications
  and presentation

In your project which uses Seis.jl, you should add then import one or more of
these as shown above.

The figures produced here were made using CairoMakie unless otherwise noted.

!!! note
    In the documentation below, wherever you see `Makie`, this represents
    whichever Makie backend module you have loaded.  In some cases, Seis.jl
    extends Makie functions with methods for Seis.jl objects

## Plotting functions

### `Makie.plot`: Default plot
Seis.jl extends `Makie.plot` to produce the default kind of plot for
whatever Seis.jl types are passed to it.  At present, it will forward to
`plot_traces` for single or array of `AbstractTrace`s.

### `plot_traces`: Single-trace plots
`plot_traces` creates a set of single-trace plots where the independent (x)
axis is time, and the dependent (y) axis is the value of the trace,
similar to simply plotting `trace(t)` against `times(t)`.

However, `plot_traces` when applied to a single `Trace` or an
`AbstractArray{<:Trace}` offers options to scale plots relatively,
show pick times, and so on.

#### Examples
1: A simple plot for a single trace, showing the picks.

```@setup plotting_makie
using CairoMakie
```

```@example plotting_makie
using Seis, CairoMakie
t = sample_data()
plot_traces(t)
```

2: Plot the vertical components for an earthquake with the distance
   in km and station code labelled.

```@example plotting_makie
t = sample_data(:regional)
z = filter(x -> endswith(x.sta.cha, "z"), t)
labels = z.sta.sta .* ", ∆ = " .* string.(round.(distance_km.(z), digits=1)) .* "km"
plot_traces(z; label=labels, sort=:dist)
```

3: Same again, but showing the falloff in amplitude with distance by
   setting all plots to have the same scale with the `ylims` option:

```@example plotting_makie
plot_traces(z; label=labels, sort=:dist, ylims=:all)
```

#### Full docstring
```@docs
Seis.plot_traces
```


### `section`: Record sections
`section` plots a 'record section', or a set of traces where a trace's
position on the y-axis is determined by some other information.

Typically record sections show traces against distance, but `section`
supports arbitrary values or functions to place the traces.

To align traces, pass a set of values to the `align` keyword argument.

#### Examples
##### 1. Simple distance record section
Using the radial components of a regional earthquake.

```@example plotting_makie
t = sample_data(:regional)
e, n, z = t[1:3:end], t[2:3:end], t[3:3:end]
rotate_to_gcp!.(e, n)
r, t = e, n
plot_section(r)
```

##### 2. Normalised section
Same again, but with traces [`normalise`](@ref)d so that the amplitudes
are more similar.

```@example plotting_makie
plot_section(r .|> normalise)
```

##### 3. Record section aligned on a set of picks
Let us use the picks stored with the key `:A` to align our traces. 
Here we can manually set the plot limits by passing arguments via `axis`,
which are then passed to
[`Makie.Axis`](https://docs.makie.org/stable/reference/blocks/axis).

```@example plotting_makie
t = sample_data(:array)
plot_section(t; align=:A, axis=(limits=(-10, 20, nothing, nothing),))
```

##### 4. Section, sorted by distance
We can plot traces sorted by epicentral distance, but equally spaced along
the y-axis, where traces are
[`normalize`](@ref)d to show how well the peaks are aligned, with the trace
amplitudes scaled down by half (using the `zoom` keyword), and
Makie keywords to set the y-axis label.  This time we
cut the traces relative to the `:A` pick rather than setting plot limits.

```@example plotting_makie
ts = sort(t; by=distance_deg)
plot_section(cut.(ts, :A, -10, 10) .|> normalize, "index";
    align=:A, zoom=0.5, axis=(ylabel="Trace index, increasing distance",))
```

##### 5. Section aligned on predicted phase arrival time
If you have installed [SeisTau.jl](https://github.com/anowacki/SeisTau.jl),
then adding predicted travel times and aligning on these becomes easy.

Here we calculate predicted times for the PKIKP arrival using `SeisTau`,
called via [`add_picks!`](@ref).

```@example plotting_makie
using SeisTau
add_picks!.(t, "PKIKP")
plot_section(cut.(t, :PKIKP, -10, 10); align=:PKIKP)
```

#### Interactions
If you are using an interactive Makie backend like GLMakie, when the
plot is active you can use the standard Makie interactions to change the plot:

- scroll to zoom
- click and drag to set a new view region
- right-click and drag to move the viewport
- hold down the `x` or `y` key while scrolling or dragging to limit the
  change to be just horizontal or vertical, respectively
- Ctrl-click to return the view to the original plot zoom and limits

In addition, with `plot_section` you can press:
- `=` to increase the trace amplitude
- `-` to decrease the trace amplitude

#### Downsampling
When scrolling and zooming with a `plot_section` plot, you will also note
that if there are too many data points to plot, the small text 'Downsampled'
appears at the bottom right of the figure.  `plot_section` automatically reloads
data into the window at either full resolution or downsampled to keep
interactions quick.  Use the keyword argument `decimate=false` to turn off
downsampling.

!!! warning
    If using `decimate=false` with very large datasets (many millions of
    points) plotting may become extremely slow or even unstable.

#### Full docstring
```@docs
Seis.plot_section
```
