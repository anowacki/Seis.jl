# Seis.jl specification

An open, fast and flexible framework for analysing seismic data in Julia.

That is the objective of Seis.jl, a future package for Julia which enables
processing of continuous passive seismic data.

If you want to help create Seis.jl, please
[get in touch](mailto:a.nowacki@leeds.ac.uk).

## Motivation

[Julia](https://julialang.org) is perhaps the most promising numerical
programming language available to the scientific community today, offering
C-like speed, inbuilt multi-dimensional (complex) arrays, parallel processing
and other goodies, all in a dynamic environment that makes writing software
quicker and less error-prone.

Passive seismologists would like to be able to take advantage of Julia to
process seismic data, but at present there isn't a widely-used package to
which we can all contribute, something which is necessary to avoid duplication
of effort.

There are a couple of reasons that a community Julia package for seismic
processing has yet to emerge:

1. Time series analysis is in many ways straightforward, so implementing
   the algorithms you need on an ad-hoc basis isn't too hard.  This way you
   have complete control individually on your workflow.
1. There are already a lot of other processing and analysis packages, so if
   you want to do something already implemented, then why bother writing
   something new?
1. Julia is a young language and is only soon to have its 1.0 release.

Despite this, there are some features of a passive seismic processing workflow
that would be desirable but aren't yet widely available:

1. Ease of extending existing functionality.  Current packages are not easy for
   newcomers to adapt to their needs for a few reasons.  Established codes
   like SAC and Seismic Handler are implemented in static compiled languages
   which though fast to run, are slower to develop in.  Newer codes like
   ObsPy are written in dynamic languages which, though fast to develop in,
   are slow to run when written natively; instead core functionality must be
   written in C or Fortran to be quick.
1. The underlying software (i.e., the language or runtime) should be
   open-source so that scientific results are reproducible.  This rules out
   implementations based on closed-source languages such as Matlab.

We believe Julia has the capability to address these needs because:

1. Julia code can run as fast (or faster!) than C or Fortran, yet is dynamic.
   In our experience, though first-pass implementations of algorithms are often
   quick anyway in Julia, further speed-up to Fortran-like speed is
   always possible within the language, without having to rewrite the code in
   something like C or Fortran.
1. Julia is open-source and maintained by both an active community of users,
   as well as the original creators.

## Scope

Seis.jl is intended to facilitate passive seismology.  That is, it intends
to provide ways to use seismic data acquired by passive seismic instruments
first and foremost.  Its current paradigm that each seismometer is a kind of
point measurement recording continuously in time.  Data are typically analysed
to understand Earth structure, tectonic processes, volcanic activity, seismicity
and many other problems.  Although active-source seismic data recorded on
geophone arrays can be processed with Seis.jl, this is not the primary focus.
Similar projects to Seis.jl are SAC, Seismic Handler and ObsPy.

## Previous work

### Passive seismic packages

These directly address the situation we will implement and are learning
exercises for how to design Seis.jl.

- [Andy Nowacki](http://homepages.see.leeds.ac.uk/~earanow)'s
  [SAC.jl](https://github.com/anowacki/SAC.jl) package offers a processing
  workflow based around SAC.  This implements much of SAC's basic signal
  processing functionality, using other Julia packages, such as filtering,
  rotation, great circle computation.  Potting is provided by
  [SACPlot.jl](https://github.com/anowacki/SACPlot.jl).
- [Joshua Jones](https://github.com/jpjones76)'s
  [SeisIO](https://github.com/jpjones76/SeisIO.jl) package offers a generic
  geophysical time series data workflow.  Basic merging, sorting and extraction
  are implemented.  Clients for FDSN and IRIS data are implemented.  Reads
  miniSEED, SAC, SEG-Y, P-SEG-Y, UW and Win32 formats.

### Previous work on active-source/reflection seismics

These are not facing passive/continuous data, but may inspire some of our
thinking in designing Seis.jl.

- [U Alberta's Signal Analysis & Imaging Group](http://saig.physics.ualberta.ca)
  offer a range of tools for reflection seismology as part of their
  [Seismic Julia](https://github.com/SeismicJulia) GitHub organisation.
- UBC's [Slim Group](https://www.slim.eos.ubc.ca) offer some codes
  which appear aimed at waveform inversion of reflection seismic data
  (e.g., [their WAVEFORM package](https://github.com/slimgroup/WAVEFORM.jl)).
- Others?

## Specification

This details what Seis.jl should do.  An example of the sort of specification
we may arrive at is:

- Provide types for traces, which includes information about the station.
  This should include:
  - Station information:
    - Name, network, location id
    - Geographic location (lon, lat, altitude, burial depth)
    - Sampling interval
    - Date and time of samples
  - Component information:
    - Azimuth, inclination, name
	- Response?
  - Anything else (metadata; implemented as `Dict` or similar)
- Provide types for events (or place this in the trace type)?
  - Lon, lat, depth, origin time, id
  - Moment tensor?
  - Anything else (metadata)
- Implement functionality (possibly in secondary packages?)
  - IO.  At least formats already implemented in other packages
    (SAC, SEG-Y, SEED)
  - Standard DSP (filtering, instrument response removal, etc.)
  - Trace manipulation (rotation, trend removal, etc.)
  - Plotting: individual and multiple traces, hodograms, record sections...
  - Direct download of data from IRIS and so on.
  - Travel times and paths in 1D Earth model
  - Array analysis (stacking, vespagrams, FK spectra, etc.)

The base Seis.jl package should probably not:

- Implement things like event location, velocity inversions, and
  anything beyond basic IO and processing.  This is what a seismic observatory
  package (SeisObs.jl?) might want to do.
- Rely too heavily on large external dependencies.  We want Seis.jl to be
  lightweight enough to use everywhere.

### Open questions

- How much do we prescribe within the trace structure?  I.e., should we provide
  fields for as many common data as we think are needed (e.g., moment tensor)
  so that these are always used, or allow a consensus to emerge on metadata
  fields?  Do we fear everyone using a `meta` field to store the same information
  under different names?
- Do we provide types for collections of traces associated with one event,
  or just operate on arrays of traces, each with own event information?  The
  latter is more akin to SAC, the former ObsPy.
- Only allow continuous traces, or allow gaps?
- How abstract should the trace representation be?  One can imagine a trace
  simply providing a lazy view into a database-like store (good for continuous data
  and e.g. ambient noise processing or event scanning), or we can just use `Array`s.  
- Do we let users choose precision and use own Array types (i.e., parameterise
  the structures on these like `struct Trace{T<:Number,V<:AbstractVector{T}}`) or
  prescribe them?

### Other implementation thoughts

- Plotting is still Julia's main weakness.  Plots.jl is probably the best way
  to go for now, but there is still no consensus plotting package.  Makie.jl
  should be compatible with Plots.jl recipes.
- We should implement e.g. filtering and response removal natively immediately.
- But rely on existing code initially for more complicated problems, e.g. travel
  time calculation.  Look into JavaCall for TauP (though I know nothing about
  Java).

## Wider questions

- How do we keep the managing of Seis.jl sustainable whilst actually doing
  research/teaching?
- How do we foster contributions from anyone using Seis.jl around the world?
  Seis.jl should be useful to all global/regional seismologists, so will begin
  with gaps where we don't do certain things.
- How do we decide what goes in Seis.jl, and what should be elsewhere?  I.e.,
  how do we review the scope?
- Which licence is appropriate?

## Members

Everyone is welcome to participate in the development of SAC.jl.  If you are
interested, please [contact Andy Nowacki](mailto:a.nowacki@leeds.ac.uk).

At this stage, we are formulating the ideas which will ensure Seis.jl is useful
to the whole passive seismic community going forward.  Then we will write it
and make it available under an open-source licence.
