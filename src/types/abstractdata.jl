"""
    AbstractData

Abstract supertype of all single-channel datatypes in Seis.jl.  An
`AbstractData` represents some data acquired at a point where it
makes sense to associate that recording with a channel, and where the
recording is limited to some period of time.  However, the recording
may be in the time or frequency domain, and need not be stationary or
evenly-sampled.
"""
abstract type AbstractData end
