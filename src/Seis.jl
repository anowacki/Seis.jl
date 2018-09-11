"""
# `Seis` – Passive seismology in Julia

## Types

Seis exports three important types:

- `Event`, which holds information about a seismic event;
- `Station`, holding information about a seismic station or channel; and
- `Trace`, which combines the above two with a time series.

In each case, the user is allowed to access some of the fields, because
`getproperty` and `setproperty!` are defined for arrays of these types.
Operating on arrays of `Trace`s is the fundamental design of Seis and offers
the greatest flexibility in dealing with seismic data, including using Julia's
broadcasting on these types with the methods defined here.
"""
module Seis

export
    # Types
    Station,
    Event,
    Trace,
    # 'Getters'
    dates,
    nsamples,
    picks,
    times,
    trace,
    # 'Setters'
    add_pick!,
    add_picks!,
    clear_picks!,
    # Geometry
    azimuth,
    backazimuth,
    distance_deg,
    distance_km,
    # Travel times
    travel_time,
    # IO
    SACtr,
    read_sac,
    write_sac

import Missings
import Missings: Missing, ismissing, missing
using Compat.Dates
import MacroTools: @capture

import SAC
import SAC: SACtr
import Geodesics
import TauPy

include("types.jl")
include("show.jl")
include("input_output.jl")
include("geometry.jl")
include("traveltimes.jl")
include("util.jl")

end # module
