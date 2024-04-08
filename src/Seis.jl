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

- `Seis.Plot`: Plotting routines.
- `Seis.Synth`: Routines to create simple wavelets.
"""
module Seis

export
    # Types
    AbstractFourierTrace,
    AbstractTrace,
    CartEvent,
    CartStation,
    CartTrace,
    Event,
    FourierTrace,
    GeogEvent,
    GeogStation,
    Station,
    Trace,
    # 'Getters'
    are_orthogonal,
    dates,
    enddate,
    endtime,
    frequencies,
    is_east,
    is_north,
    is_horizontal,
    is_vertical,
    nearest_sample,
    nfrequencies,
    nsamples,
    picks,
    startdate,
    starttime,
    times,
    trace,
    # 'Setters'
    add_pick!,
    add_picks!,
    clear_picks!,
    origin_time!,
    origin_time,
    # Geometry
    azimuth,
    backazimuth,
    distance_deg,
    distance_direct,
    distance_km,
    incidence,
    # Operations
    cut!,
    cut,
    decimate!,
    decimate,
    differentiate!,
    differentiate,
    envelope!,
    envelope,
    fft,
    flip!,
    flip,
    ifft,
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
    resample!,
    resample,
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
    read_mseed,
    read_sac,
    write_mseed,
    write_sac,
    write_sac_header,
    # Sample data
    sample_data,
    # Rotation
    rotate_through!,
    rotate_through,
    rotate_to_azimuth_incidence!,
    rotate_to_azimuth_incidence,
    rotate_to_enz!,
    rotate_to_enz,
    rotate_to_gcp!,
    rotate_to_gcp,
    rotate_to_lqt!,
    rotate_to_lqt,
    sort_traces_right_handed,
    # Plotting
    plot_traces,
    # Analysis
    spectrogram

using Dates
using LinearAlgebra
using LinearAlgebra: ×, ⋅, norm
using Statistics: mean, covm, varm

import Glob
import DSP
import DSP: resample, spectrogram
import FFTW
import FFTW: fft, ifft
import StaticArrays

import Geodesics

include("compat.jl")

# All basic types
include("types.jl")
include("spec.jl")
include("conversion.jl")

# File formats submodules
include("io/SAC/SAC.jl")
using .SAC
include("io/Miniseed/Miniseed.jl")
using .Miniseed

include("show.jl")
include("input_output.jl")
include("geometry.jl")
include("traveltimes.jl")
include("util.jl")
include("operations.jl")
include("rotation.jl")
include("filtering.jl")
include("sample_data.jl")
include("spectrogram.jl")

# Functionality submodules
include("Synth.jl")
import .Synth

# Plotting functionality is only loaded when using Plots or Makie
include("Plot.jl")
include("makie.jl")

end # module
