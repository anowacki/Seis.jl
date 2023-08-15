# Seis.jl v0.4 release notes

## Julia compatibility
- Seis.jl now supports Julia v1.6 (long-term support release) upwards.
  Versions older than this are not supported.

## Breaking changes
- `rotate_through[!]` now explicitly implements the partial de-facto
  previous behaviour: that the direction of the rotation is defined
  by the order of the traces.
  **Versions of this function prior to this change may have got the
  sense of the rotation wrong and defined the trace azimuths wrong**
  depending on the trace input order,
  and so this is technically a breaking change.  However the new
  behaviour is correct, so can now be relied upon.
- `read_mseed`:
  - no longer supports the `maximum_offset` keyword argument;
  - moves the positional argument `T` specifying the trace type to the
    end to better match up with `Base.read`; and
  - now returns a `Trace{Float64, Vector{Float32}, Seis.Geographic{Float64}}`
    by default.
- `plot`: The deprecated `picks` keyword argument has been removed in
  favour of `show_picks`.

## New features and non-breaking changes
### New types
- There is a new `FourierTrace` type, which contains trace data
  in the frequency domain.  Most operations work on this type as for
  time-series `Trace` data.  Get a `FourierTrace` from a `Trace` by calling
  `fft` on it, and go back with `ifft`.
### IO
- You can now write miniSEED files with `write_mseed`.
- `read_sac` and `write_sac` will accept an `IO` object (like an
  `IOBuffer`) to read and write to.
- The `echo` keyword argument to `read_sac(glob, dir)` now defaults to
  `false`, meaning matching files are no longer written to stdout
  unless specifically requested.
- `write_sac_header` will write just the header part of a `Trace` to a
  SAC file.  This is useful when reading just the headers to update
  header information on disk without reading and writing the whole
  trace data.
### Trace operations
- `resample[!]` can resample traces to arbitrary sampling intervals.
- `rotate_through[!]` can rotate arbitrary pairs of orthogonal traces,
  not just horizontal ones.
- `rotate_to_enz[!]` rotates triplets of orthogonal components to ENZ
  orientation
- `rotate_to_lqt[!]` rotates triplets of orthogonal components to
  LQT orientation.
- `rotate_to_azimuth_incidence[!]` rotates triplets of orthogonal components
  to arbitrary orientations.
- Spectrograms can be calculated with `spectrogram`, and can be plotted
  with `Seis.Plot.plot_spectrogram`.
### Picks
- `add_picks!` now accepts an absolute time (a `DateTime`) at which to
  record a pick
### Convenience functions
- The new single-argument `origin_time` method returns the origin time of
  the trace (like `.evt.time` but nicer-looking).

## Deprecated or removed
- `traces_are_orthogonal` will be removed in v0.5.  It has been replaced by
  `are_orthogonal` which works for non-horizontal pairs as well.

## Notable bug fixes
- `rotate_through[!]` behaved in an inconsistent way and has been fixed
  (see 'Breaking changes'), but its new behaviour is different.
