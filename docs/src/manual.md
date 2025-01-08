# User manual

This section describes how to use Seis to perform various basic
processing operations on seismic data.

## Preamble

The examples in this manual all assume that you have first used
the module like so:

```@repl example
using Seis
```

## Sample data
Seis comes with a small selection of sample data which will be used for
the examples in this manual.  See the help for [`sample_data`](@ref)
for details.

## Types
Working with Seis requires understanding a little bit about the types
it uses.  Fundamentally, there are four important types:
- [`Trace`](@ref)s, which contain
- [`Station`](@ref)s and
- [`Event`](@ref)s; and also
- [`FourierTrace`](@ref)s, which are frequency-domain `Trace`s.

Each of these types have fields which either must contain data,
or allow data to be [`missing`](https://docs.julialang.org/en/v1/manual/missing/#missing-1).

### `Trace`s
A `Trace` holds a single continuous record of evenly-sampled data.
As such, it cannot contain gaps or overlaps.  It may be referenced
to an absolute start time, and by convention all times of day in Seis
are [UTC](https://en.wikipedia.org/wiki/Coordinated_Universal_Time).
Equally, however, a trace may have no absolute reference time.

Traces have a small number of fields which are directly accessible
by users:
- `delta` is the samping interval (1 / sampling rate) in s;
- `b` is the start time of the trace, relative to any absolute time
  if any;
- `sta` is a [`Station`](@ref Stations) holding information about the
  site where this recording was made;
- `evt` is an [`Event`](@ref Events) holding information about any
  event related to this recording, such as an earthquake;
- `picks` holds a set of times of importance such as arrival picks; and
- `meta` holds any other information.

None of these fields can be `missing`.

To access the data stored in a trace, use the [`trace`](@ref) function
to return the underlying array.

The user-accessible fields can be retrieved or modified as usual for
composite types:

```@repl example
# A trace starting a 0 s with sampling interval 1 s made of 100 random samples
t = Trace(0, 1, rand(100))
# Update the sampling interval to 0.1 s
t.delta = 0.1
```

#### Mutability

`Trace`s are mutable.  This means that the values they hold can be
changed by the user at will.  It is often useful to copy a trace after
you have performed an operation on it.  Thoughout Seis, as with Julia
generally, functions which modify a `Trace` have an exclamation mark
(`!`) at the end of their name.

!!! note
    `Trace`s are 'parameterised' on the type of floating point numbers
    they use, the type of data storage, and the geometry in which they
    are defined (i.e., either in geographic or Cartesian coordinates).
    A `Trace` is technically a `Trace{T, V, P} where {T, V, P}`, where
    `T` is the number type used, `V` is the type of vector holding the
    trace data, and `P` is the geometry type.
    This is why a trace has a scary-looking type like
    `Trace{Float32, Array{Float32, 1}, Seis.Geographic{Float32}}`.
    You only need to deal with these type parameters if you want to
    change these defaults.  See [Geometry](@ref) for more information.

### Collections of `Trace`s
In Seis, methods are usually written to accept `Trace`s.  There is no
special container for several traces; instead you can use `Array`s of
traces just as you would any other collection of types.

For example, to find the epicentral distance between the earthquake
and a number of recordings of it at different stations, you could
use Julia's broadcasting feature to call the [`distance_deg`](@ref)
function on each `Trace`.  (Note that the `sample_data` function
returns a simple `Vector{Trace}` of recordings.)

```@repl example
t = sample_data(:array)
typeof(t)
eltype(t)
distance_deg.(t)
```

To find the nearest station, you could write the following:

```@repl example
t[argmin(distance_deg.(t))].sta
```


### `Station`s
`Station`s define a stationary point in space where a recording was made.
They have the following fields which you can access and modify:

- `net`, `sta`, `loc` and `cha` are used respectively for the network,
  station, channel and location names of a single recording channel.
  They correspond to the
  [SEED channel naming convention](http://www.fdsn.org/pdf/SEEDManual_V2.4_Appendix-A.pdf), but all can be `missing`.
- `lon`, `lat` and `elev` are the longitude and latitude (in °, positive
  eastwards and northwards) of the station, whilst `elev` is the elevation
  above sea level in m.
- `azi` and `inc` together define the orientation of the channel; that is,
  in which spatial direction positive values point.  `azi` is the local
  azimuth in ° measured clockwise from north, whilst `inc` is the local
  inclination in ° measured from the vertical downwards.  For seismic data,
  vertical channels will have `inc == 0` and `azi == 0`, whilst horizontal
  channels will have `inc == 90` and `azi == 0` for north and `azi == 90`
  for east.
- `meta` holds all other information.

See [Geometry](@ref) for more details on other coordinate systems in which
objects can be placed in Seis.

### `Event`s

An `Event` marks the source of a `Trace`.  Typically, it represents the
source of energy for a recording, or it may simply indicate the start
time of a section of data with no energy source implied.

`Event`s have the following accessible fields, all of which can be `missing`:

- `time` is a `Dates.DateTime` giving an absolute date and time in UTC
  against which a `Trace`'s `b` field (beginning time) is relative.
  For an earthquake or other source of energy, it should be the origin
  time of that event.
- `lon`, `lat` and `dep` are the longitude and latitude (in °) and
  depth (in km) of any identified source for the data, such as an earthquake.
- `id` is a `String` giving some identifying information about an
  event, and could be a catalogue ID or otherwise.
- `meta` holds all other information.  By convention, the following
  fields in `meta` might be used:
  - `mag` holds an event magnitude;
  - `quakeml` holds information about an event in QuakeML form, if any;
  - `catalog` gives the name of the catalogue for this event.

### `FourierTrace`s

A `FourierTrace` is a frequency-domain version of its corresponding `Trace`.
Usually you convert an existing `Trace` to the time domain by calling
[`fft`](@ref fft(::Trace)), and then call [`ifft`](@ref ifft(::FourierTrace))
to convert it back.  Because `Trace`s have to have real-valued data, a
`FourierTrace`'s data is a one-sided Fourier transform.  You access the underlying
data in the same way as for `Traces`, with [`trace`](@ref).

A `FourierTrace`'s public fields are:

- `delta`, giving the frequency spacing of the data;
- `evt`, giving the event associated with the recording;
- `sta`, giving station information; and
- `meta` holding other things.

`evt`, `sta` and `meta` are brought across from a `FourierTrace`'s
originating `Trace` when constructed in the usual way with `fft`,
and they are passed back to the new trace when calling `ifft`.

`FourierTrace`s have the following accessor methods of their own:

- [`frequencies`](@ref) gives the frequencies of each data point;
- [`nfrequencies`](@ref) gives the number of data points; and
- [`nsamples`](@ref) gives the number of data points in the equivalent
  time-domain `Trace`.

Because a `FourierTrace` retains enough information about the starting
trace to turn it back, you can also call [`starttime`](@ref) on a `FourierTrace`.


## Getting data in

Data can be loaded into Julia via a number of means:

- Reading data from disk in SAC or miniSEED format.

  For this use either [`read_sac`](@ref) or [`read_mseed`](@ref).

  To read a single file, use the one-argument form:

  ```
  julia> t1 = read_sac("single_file.sac")
  Seis.Trace{Float32,Array{Float32,1},Seis.Geographic{Float32}}:
             b: 52.66
         delta: 0.01
  Station{Float32,Seis.Geographic{Float32}}:
       sta.lon: -120.0
       sta.lat: 48.0
       sta.sta: CDV
       sta.azi: 0.0
       sta.inc: 0.0
      sta.meta: Seis.SeisDict{Symbol,Any}()
  Event{Float32,Seis.Geographic{Float32}}:
       evt.lon: -125.0
       evt.lat: 48.0
       evt.dep: 0.0
      evt.time: 1981-03-29T10:38:14
        evt.id: K8108838
      evt.meta: Seis.SeisDict{Symbol,Any}()
  Trace:
         picks: 2
          meta: SAC_lpspol => true
                SAC_nevid => 0
                SAC_iftype => 1
                file => "../test/test_data/seis.sac"
                SAC_idep => 50
                SAC_iztype => 9
                SAC_lcalda => true
                SAC_unused18 => false
                SAC_lovrok => true
                SAC_norid => 0
                SAC_ievtyp => 42
  
  julia> t2 = read_mseed("file.mseed")
  2-element Array{Trace{Float32,Array{Float32,1},Seis.Geographic{Float32}},1}:
   Seis.Trace(GB.CWF..BHZ: delta=0.02, b=0.0, nsamples=3000)
   Seis.Trace(GB.CWF..HHZ: delta=0.01, b=0.0, nsamples=6000)
  ```
  
  !!! note
      SAC files can contain only one single continuous data channel, whilst
      miniSEED files can contain more than one, and so an array of `Trace`s
      are returned.

  Reading several files using a globbing pattern:

  ```
  julia> read_sac("Event_*/*Z.sac", "DATA")
  ```

- Downloading data from a remote server.

  For this, install and use
  [SeisRequests](https://github.com/anowacki/SeisRequests.jl).

  ```
  julia> using SeisRequests

  julia> t = get_data(code="IU.ANMO.00.BH?", starttime="2018-02-02", endtime="2018-02-02T01:00:00") # an hour of data
  [ Info: Request status: Successful request, results follow
  3-element Array{Trace{Float64,Array{Float64,1},Seis.Geographic{Float64}},1}:
   Seis.Trace(IU.ANMO.00.BH1: delta=0.05, b=0.0, nsamples=72000)
   Seis.Trace(IU.ANMO.00.BH2: delta=0.05, b=0.0, nsamples=72000)
   Seis.Trace(IU.ANMO.00.BHZ: delta=0.05, b=0.0, nsamples=72000)
  ```


## Creating data from scratch

[`Trace`](@ref)s can be constructed using the constructors.  For example:

```@repl example
b = 0
delta = 0.01
data = randn(1000)
t = Trace(b, delta, data) # Fill trace with data already available
```


## Writing data out

Seis currently supports writing in SAC and miniSEED format,
using respectively the [`write_sac`](@ref) and [`write_mseed`](@ref) functions.

```@repl example
t = sample_data()
write_sac(t, "outfile.sac")
write_mseed("outfile.mseed", t)
```

Note that miniSEED files can contain several traces, and so you can also
pass an array of traces to `write_mseed`:

```@repl example
t = sample_data(:local)
write_mseed("local_data.mseed", t)
```

!!! note
    `Trace`s written in miniSEED format must have their `.evt.time` field
    set, since this is required in the file format.

!!! note
    miniSEED files do not contain any information about the station or
    event.  The only information saved is the station network, station,
    location and channel codes and the start date of the first sample.

!!! note
    The order of arguments (file name and trace object) is different
    between `write_mseed` and `write_sac` for legacy reasons.  In a
    future major update the functions will all have the same order as
    `Base.write` (i.e., file name or `IO` object first, then trace
    object).

Seis supports reading and writing only the header part of SAC files:
- To read just headers, use `read_sac(file; header_only=true)`.
- To write just headers, use [`write_sac_header`](@ref).

miniSEED files can be read with only headers using
`read_mseed(file; header_only=true)`.


## Basic processing

Once read in or obtained some other way, data can be processed in a number
of ways.

Note that, as mentioned above ([Mutability](@ref)), traces can be modified
in place, or copies taken when processing steps are applied.  Modifying in
place is usually faster as it does not involve copying all the information,
but it is often more convenient to copy traces as well.  For this reason,
there is always a modifying version of a function (ending in `!`), and a
copying version, without.

- Cut the data to a certain window, in either relative time, absolute
  time, or relative to a time pick, with [`cut`](@ref) or [`cut!`](@ref).
- Decimate the data, with or without an antialiasing filter
  ([`decimate`](@ref), [`decimate!`](@ref)).
- Remove the mean value of the trace with [`remove_mean`](@ref) or
  [`remove_mean!`](@ref).
- Remove a linear trend in the data with [`remove_trend`](@ref) or
  [`remove_trend!`](@ref).
- Taper the signal with [`taper`](@ref) or [`taper!`](@ref).
- Filter the data, performing a low-pass ([`lowpass`](@ref),
  [`lowpass!`](@ref)), high-pass ([`highpass`](@ref), [`highpass!`](@ref)),
  band-pass ([`bandpass`](@ref), [`bandpass!`](@ref)) or
  band-stop ([`bandstop`](@ref), [`bandstop!`](@ref)) filter.  For each,
  an acausal version of the filter can be obtained by passing the
  `twopass=true` option.
- Differentiate ([`differentiate`](@ref), [`differentiate!`](@ref)) or
  integrate ([`integrate`](@ref), [`integrate!`](@ref)) the trace.
- Normalise the data to a certain value with [`normalize`](@ref) or
  [`normalize!`](@ref) (or [`normalise`](@ref) for us Brits).
- Take the envelope of a trace with [`envelope`](@ref) or [`envelope!`](@ref).
- Change the trace sampling interval with [`resample`](@ref) or
  [`resample!`](@ref).
- Merge multiple traces with the same channel code together with [`merge!`](@ref)
  and [`merge`](@ref).

In addition, sets of traces can be rotated:
- [`rotate_through`](@ref) and [`rotate_through!`](@ref) rotate
  a pair of traces by a certain angle.
- [`rotate_to_gcp`](@ref) and [`rotate_to_gcp!`](@ref) rotate pairs
  to radial (pointing away from the source at the station)
  and transverse (90° clockwise from this) components.
- [`rotate_to_lqt`](@ref) and [`rotate_to_lqt!`](@ref) rotate triplets
  of orthogonal compontents to longitudinal (L), transverse (T) and
  Q (orthogonal to the others) orientations.
- [`rotate_to_azimuth_incidence`](@ref) and
  [`rotate_to_azimuth_incidence!`](@ref) rotate triplets of orthogonal
  traces to arbitrary orientations.

All of the above will use trace header information (coordinates of the
event and station, and channel orientations) to automatically compute
the directions if possible.


## Accessing raw trace data

Beyond these basic functions, `Trace`s are designed to be used to build
your own processing workflows on top.  For example, shear wave splitting
analysis can be done with [SeisSplit](https://github.com/anowacki/SeisSplit.jl)
and array processing with
[Beamforming](https://github.com/anowacki/Beamforming.jl).

To access the raw trace, use the [`trace`](@ref) function, which returns
the data:

```@repl example
t = Trace(0, 0.1, [1, 2, 3, 4, 5])
data = trace(t)
```

As is usual with Julia, `data` is a variable bound to the same object as
the underlying data of the trace `t`.  Modifying `data` will also modify the
values in `t`:

```@repl example
data[1] = 0 # Set the first point to be zero
trace(t) # Note the first point is now 0
```

`data` can be
[`empty!`](https://docs.julialang.org/en/v1/base/collections/#Base.empty!)d, or
[`push!`](https://docs.julialang.org/en/v1/base/collections/#Base.push!)ed
or
[`append!`](https://docs.julialang.org/en/v1/base/collections/#Base.append!)ed
to, and so on—it is simply a vector
of data points.

## Accessing other properties
### Times
The first sample of a trace `t` occurs at time `t.b`, and each sample is
`t.delta` seconds after the previous one.  Therefore, each sample of the
trace has a time, and this array (actually a subtype of `AbstractRange`)
can be obtained with [`times`](@ref):

```@repl example
times(t)
```

For consistency, the start and end times of a trace are given by
[`starttime`](@ref) and [`endtime`](@ref):

```@repl example
starttime(t), endtime(t)
```

All these times are relative to the event date if that is set.  If it
is, then [`startdate`](@ref) and [`enddate`](@ref) give the start and
end date of the trace in UTC:

```@repl example
using Dates
t.evt.time = DateTime(2013, 5, 24, 5, 44, 49, 880) # Set the origin time
startdate(t), enddate(t)
```

Like before, we can get a set of dates for every sample with [`dates`](@ref):

```@repl example
dates(t)
```

There should be 5 values here, since we gave the trace a set of 5 data points.
[`nsamples`](@ref) tells us this:

```@repl example
nsamples(t)
```

The closest sample to a certain time or date can be found with
[`nearest_sample`](@ref):

```@repl example
nearest_sample(t, 0.11)
nearest_sample(t, DateTime(2013, 5, 24, 5, 44, 50))
```


## Picks
Seis allows you to assign arbitrary time picks to a `Trace` and use
these in flexible ways.

Picks are simply markers which label a certain point in time with a label.
They are useful for marking particular phase arrivals and for other
purposes.  All picks in Seis are relative to the origin time for the
trace if any.  So if a trace's `starttime` is 3 s and a pick is at 4 s,
the pick points to a time 1 s after the trace starts.

Picks have two fields: `time` and `name`.  `time` is the time in seconds
after zero time, which is the event time if set.  `name` may be a `String`
such as `"PcP"`, or `missing` if no name is needed for this pick.

### Getting picks
Picks are accessible in the `picks` field of a trace:

```@repl example
t = sample_data();
t.picks
```

You can see that there are two picks defined: one at 53.67 s with key `:A`,
and one at 60.98 s with key `:F`.  Each pick is associated with a key, which
can be a `Symbol`, or an `Int`.  Picks with a `Symbol` key are named picks and
can be accessed in two ways.

You can access a named pick with dot notation if you know the literal key:

```@repl example
t.picks.A
```

You can also use `getindex` or `[]` notation to access a named pick using
a variable bound to a `Symbol` or a literal value:

```@repl example
t.picks[:A]
pickkey = :A; t.picks[pickkey]
```

For any pick, access its time and name like so:

```@repl example
t.picks.A.time, t.picks.A.name
```

### Setting picks
Setting picks is done similarly.  To add a named pick with name `:S`,
time 2 s and name `"Sg"`, you can do

```@repl example
t.picks.S = 2, "Sg"
t.picks
```

Or, you can use `setindex!` (`[]`):

```@repl example
t.picks[:S] = 2.5
t.picks
```

Note that the `:S` pick is overwritten.  Note also here that by not providing
a name for the pick, its name defaults to `missing`.

Picks can also be set using [`add_picks!`](@ref), which can take a
`Dates.DateTime` to set a pick in absolute time rather than a time relative
to the trace origin.

```@repl example
using Dates
add_pick!(t, DateTime("1981-03-29T10:39:10"), "Coffee time")
```

### Removing picks
Remove picks as usual for `Dict`-like collections of things, with `delete!`:

```@repl example
delete!(t.picks, :S)
t.picks
```

You can also use the special method for the `picks` field of setting a
pick to `missing`:

```@repl example
t.picks.F = missing
t.picks
```

Remove all picks with the [`clear_picks!`](@ref) function:

```@repl example
clear_picks!(t)
t.picks
```

### Arrays of picks
Because `t.picks` is a special type of `Dict` (a [`Seis.SeisDict`](@ref)),
we can also access the `time` and `name` fields of sets of picks easily.

For example, if dealing with a set of stations where we have a pick
for each, we can extract the vector of all pick times simply:

```@repl example
t = sample_data(:array)
t[1].picks # Picks of the first trace
t.picks.A.time # Times of all A picks for all traces
```



## Geometry
Seis allows you to have traces defined in any coordinate system you like.
By default, it provides two main types:
- `Geographic`, where coordinates are given as longitudes and latitudes
  in °, and depths in km below sea level, and
- `Cartesian`, where coordinates are ``x``, ``y`` and ``z`` in a right-handed
  system and each is given in m.  Conventionally in Seis, ``x`` is a local
  easting and ``y`` is a local northing, meaning ``z`` is the upwards direction
  and usually ``z`` is ``0`` at sea level or the local reference level.

By default in Seis, all `Event`s, `Station`s and `Trace`s are geographic,
and use the `Geographic` type.  Therefore you do not need to do anything
special to work in geographic sense with longitudes, latitudes, and so on.

### Cartesian geometry
However, if you are dealing with data where it makes more sense to
consider positions in a Cartesian system, that is straightforward.
Use the special constructors for Cartesian objects to create them:
- [`CartTrace`](@ref) creates a new `Trace` in an ``x, y, z`` system.
- [`CartStation`](@ref) creates a Cartesian-referenced `Station`.
- [`CartEvent`](@ref) makes a Cartesian-referenced event.

Where it makes sense, certain accessor functions like [`distance_km`](@ref)
are defined for these types, and these are necessarily different to
those for geographic objects.

For example, [`distance_direct`](@ref) will calculate the straight-line
distance in m between the event and station when in Cartesian coordinates:
```@repl example
t = CartTrace(0, 1e-4, 1000)
t.sta.x, t.sta.y, t.sta.z = 101, -35, 12
t.evt.x, t.evt.y, t.evt.z = 12, 14, -100
distance_direct(t)
```

### Geometry type hierarchy and custom geometry types
`Geographic` and `Cartesian` are subtypes of [`Seis.Position`](@ref).
This means that if you want to define your own coordinate system,
you are able to do so by subtyping `Seis.Position` and you can then write
methods taking a `Trace{T, V, P} where {T, V, P}`, where `P` is your new
type.  These will replace the generic methods allowing special processing
to be done in cases where the generic calculations need to be changed,
without requiring writing duplicate methods where they do not.
An example of a different coordinate system might be when stations
have [positions which vary in time](http://geoweb.princeton.edu/people/simons/earthscopeoceans/).

### [`Seis.Geographic`](@ref)

### [`Seis.Cartesian`](@ref)
