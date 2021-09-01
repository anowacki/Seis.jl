# Function index

## Public types and functions
### Types

```@docs
Trace
AbstractTrace
CartTrace
Event
CartEvent
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
is_east
is_north
is_horizontal
is_vertical
nearest_sample
nsamples
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
origin_time
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
flip!
flip
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
rotate_to_gcp!
rotate_to_gcp
```

### IO
```@docs
read_mseed
read_sac
write_sac
Seis.SAC.SACTrace
```

### Example data sets
```@docs
sample_data
```
