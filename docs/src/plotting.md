# Plotting

## Introduction
Seis comes with plotting functionality which relies on
[Plots](https://docs.juliaplots.org/latest/).  To create plots with
Seis, first install Plots:
```
import Pkg
Pkg.add("Plots")
```

Then whenever you wish to create plots for `Trace` objects, do
```
julia> using Seis.Plot, Plots
```

This will bring the recipes into scope.

## 'Recipes'
Seis.jl implements plotting via so-called
'[recipes](https://docs.juliaplots.org/latest/recipes/)', based on the
[Plots](https://docs.juliaplots.org/latest/) plotting package.
Recipes allow us to define plots for `Trace`s without requiring the
user to install Plots.  However, if Plots is installed and a
`using Plots` command has been issued, then plots for Seis's
objects can be easily created.

## Plotting functions

### `plot`: Single-trace plots
`plot` creates a set of single-trace plots where the independent (x)
axis is time, and the dependent (y) axis is the value of the trace,
similar to simply plotting `trace(t)` against `times(t)`.

However, `plot` when applied to a single `Trace` or an
`AbstractArray{<:Trace}` offers options to scale plots relatively,
show pick times, and so on.

#### Examples
1: A simple plot for a single trace, showing the picks.

```@setup plotting
using Plots
default(fontfamily="Helvetica")
```

```@example plotting
using Seis, Seis.Plot, Plots
t = sample_data()
plot(t)
```

2: Plot the vertical components for an earthquake with the distance
   in km and station code labelled.

```@example plotting
t = sample_data(:regional)
z = filter(x -> endswith(x.sta.cha, "z"), t)
labels = z.sta.sta .* ", ∆ = " .* string.(round.(distance_km.(z), digits=1)) .* "km"
plot(z, label=labels, sort=:dist)
```

3: Same again, but showing the falloff in amplitude with distance by
   setting all plots to have the same scale with the `ylims` option:

```@example plotting
plot(z, label=labels, sort=:dist, ylims=:all)
```

#### Full docstring
```@docs
Seis.Plot.plot
```


### `section`: Record sections
`section` plots a 'record section', or a set of traces where a trace's
position on the y-axis is determined by some other information.

Typically record sections show traces against distance, but `section`
supports arbitrary values or functions to place the traces.

To align traces, pass a set of values to the `align` keyword argument.

#### Examples
1: Simple distance record section for the radial components of a
   regional earthquake.

```@example plotting
t = sample_data(:regional)
e, n, z = t[1:3:end], t[2:3:end], t[3:3:end]
rotate_to_gcp!.(e, n)
r, t = e, n
section(r)
```

2: Same again, but with traces [`normalise`](@ref)d so that the amplitudes
   are more similar.

```@example plotting
section(r .|> normalise)
```

3: Record section aligned on a set of picks stored with the key `:A`.

```@example plotting
t = sample_data(:array)
section(t, align=:A, xlim=(-10, 20))
```

4: Section, sorted by distance but equally spaced, where traces are
   [`normalize`](@ref)d to show how well the peaks are aligned, with the trace
   amplitudes scaled down by half (using the `zoom` keyword), and
   defining other Plots keywords to set the y-axis label.

```@example plotting
section(t .|> normalize, sortperm(t, by=distance_deg), align=:A,
    xlim=(-10, 10), zoom=0.5, ylabel="Trace index, increasing distance")
```

If you have installed [SeisTau.jl](https://github.com/anowacki/SeisTau.jl),
then adding predicted travel times and aligning on these becomes easy.

5: Record section aligned on the predicted PKIKP arrival time, with
   predicted times calculated using `SeisTau`, called via [`add_picks!`](@ref).

```@example plotting
using SeisTau
add_picks!.(t, "PKIKP")
section(t, align=:PKIKP, xlim=(-10, 10))
```

#### Full docstring
```@docs
Seis.Plot.section
```


### `hodogram`: Particle motion plots
Hododgrams, or particle motion plots, are parametric plots
of two components of motion through time.  These can be plotted
in Seis using `hodogram`.

#### Examples
1: Hodogram of the horizontal components of an earthquake, windowed around
   the P-wave arrival, showing the backazimuth to the event location,
   and the particle motion and backazimuth label set to have different
   colours using the `linecolor` option from `Plots`.

```@example plotting
t = e, n = sample_data(:regional)[1:2]
cut!.(t, 50, 70)
hodogram(e, n, backazimuth=true, linecolor=[:red :black])
```

#### Full docstring
```@docs
Seis.Plot.hodogram
```


### `plot_spectrogram`: Spectrogram plots
Spectrograms show the variation of frequency content with time for a
trace.  Spectrograms can be calculated with `spectrogram`, then plotted
using `plot_spectrogram`.

#### Example
```@example plotting
t = sample_data()
spec = spectrogram(t)
plot_spectrogram(spec)
```

#### Full docstring
```@docs
Seis.Plot.plot_spectrogram
```
