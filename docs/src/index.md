# Seis.jl

An open, fast and flexible framework for analysing seismic data in Julia



## What is [Seis.jl](https://github.com/anowacki/Seis.jl)?
A [Julia](http://julialang.org) package for reading and writing files
in the [QuakeML](https://quake.ethz.ch/quakeml) format, which describes
the properties of sets of seismic events, such as earthquakes and explosions.

This package is primarily meant to be used by other software to correctly
and reliably interact with QuakeML files.  For example,
[Seis.jl](https://github.com/anowacki/Seis.jl) and its related libraries
use QuakeML.jl to parse QuakeML files, but do not expose QuakeML.jl
types or functions to the user.  Though QuakeML.jl is intended to be used
as software by other software, it is still a goal that it should be easy
to use directly and well-documented and -tested.

## How to install
Seis can be added to your Julia environment like so:

```julia
julia> ] # Press ']' to enter the REPL Pkg mode
(@v1.4)> add https://github.com/anowacki/Geodesics.jl https://github.com/anowacki/Seis.jl
```

## Testing
To check that your install is working correctly, you can run the package's
tests by doing:

```julia
julia> import Pkg; Pkg.test("Seis")
```

## Contents
```@contents
```