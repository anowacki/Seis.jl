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

## Submodules

Some of Seis's functionality is not by default exported and is encapsulated into
separate modules.  To access these, do `using Seis.<submodule>`.  They are:

- `Seis.Plot`: Plotting routines which are activated when you do `using Plots`.
"""
module Seis

export
    # Types
    Station,
    Event,
    AbstractTrace,
    Trace,
    # 'Getters'
    dates,
    endtime,
    is_horizontal,
    is_vertical,
    nearest_sample,
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
    # Operations
    cut!,
    cut,
    decimate!,
    decimate,
    differentiate!,
    differentiate,
    envelope!,
    envelope,
    integrate!,
    integrate,
    normalise!,
    normalise,
    normalize!,
    normalize,
    remove_mean!,
    remove_mean,
    remove_trend!,
    remove_trend,
    taper!,
    taper,
    # Filtering
    bandstop!,
    bandstop,
    bandpass!,
    bandpass,
    highpass!,
    highpass,
    lowpass!,
    lowpass,
    # Travel times
    travel_time,
    # IO
    SACtr,
    channel_code,
    read_sac,
    write_sac,
    # Sample data
    sample_data

using Dates
using LinearAlgebra
using Statistics: mean, covm, varm

import DSP
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
include("operations.jl")
include("filtering.jl")
include("sample_data.jl")

# Submodules
include("Synth.jl")
import .Synth

# Plotting functionality is only loaded when using Plots
include("Plot.jl")

end # module
