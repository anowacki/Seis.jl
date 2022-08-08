# Function index

## Public types and functions
### Types

```@docs
Trace
AbstractTrace
CartTrace
Event
CartEvent
FourierTrace
AbstractFourierTrace
Station
CartStation
```

### Accessor functions
```@docs
are_orthogonal
channel_code
dates
enddate
endtime
frequencies
is_east
is_north
is_horizontal
is_vertical
nearest_sample
nfrequencies
nsamples
origin_time
picks
startdate
starttime
times
trace
```

### Setter functions
```@docs
add_pick!
add_picks!
clear_picks!
origin_time!
```

### Geometry functions
```@docs
azimuth
backazimuth
distance_deg
distance_direct
distance_km
incidence
```

### Trace operations
```@docs
cut!
cut
decimate!
decimate
differentiate!
differentiate
envelope!
envelope
fft(::Trace)
flip!
flip
ifft(::FourierTrace)
integrate!
integrate
normalise!
normalise
normalize!
normalize
remove_mean!
remove_mean
remove_trend!
remove_trend
resample!
resample
taper!
taper
```

### Filtering
```@docs
bandstop!
bandstop
bandpass!
bandpass
highpass!
highpass
lowpass!
lowpass
```

### Trace rotation
```@docs
rotate_through!
rotate_through
rotate_to_enz!
rotate_to_enz
rotate_to_gcp!
rotate_to_gcp
rotate_to_lqt!
rotate_to_lqt
rotate_to_azimuth_incidence!
rotate_to_azimuth_incidence
sort_traces_right_handed
```

### IO
```@docs
read_mseed
read_sac
write_mseed
write_sac
write_sac_header
Seis.SAC.SACTrace
```

### Example data sets
```@docs
sample_data
```
