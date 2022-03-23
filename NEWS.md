# Seis.jl v0.4 release notes

## Breaking changes
- `rotate_through[!]` now explicitly implements the partial de-facto
  previous behaviour: that the direction of the rotation is defined
  by the order of the traces.
  **Versions of this function prior to this change may have got the
  sense of the rotation wrong and defined the trace azimuths wrong**
  depending on the trace input order,
  and so this is technically a breaking change.  However the new
  behaviour is correct, so can now be relied upon.

## New features
- `rotate_through[!]` can rotate arbitrary pairs of orthogonal traces,
  not just horizontal ones.
- `rotate_to_lqt[!]` rotates triplets of orthogonal components to
  LQT orientation.
- `rotate_to_azimuth_incidence!` rotates triplets of orthogonal components
  to arbitrary orientations.

## Deprecated or removed
- `traces_are_orthogonal` will be removed in v0.5.  It has been replaced by
  `are_orthogonal` which works for non-horizontal pairs.


## Notable bug fixes
- `rotate_through[!]` behaved in an inconsistent way and has been fixed
  (see 'Breaking changes'), but its new behaviour is different.
