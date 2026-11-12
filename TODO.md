# TODO

List of things to do to Seis.jl

- `docs/src/manual.md:303`: Update the note about different argument
  orders between `write_mseed` and `write_sac` when this is sorted.
- Add time-limited reading for miniSEED
- Add time-limited reading for SAC including relative to headers, and
  absolute and relative time.
- Add merging of sets of traces where channels don't match, returning
  another vector of traces.  (Wait until we have a trace collection type?)