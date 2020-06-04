# [Seis.jl](https://github.com/anowacki/Seis.jl)

An open, fast and flexible framework for analysing seismic data in Julia.


## What is [Seis.jl](https://github.com/anowacki/Seis.jl)?
Seis is a [Julia](https://julialang.org) package for dealing with
seismic data.  It allows you to read and write files, and perform
processing on them.

Even more importantly, it provides a flexible foundation on which to
build your own processig workflows.  For example, the following packages
build on Seis to provide extra functionality and form the Seis ecosystem:

- [SeisRequests](https://github.com/anowacki/SeisRequests.jl) downloads
  seismic data from remote servers direct to your machine.
- [SeisSplit](https://github.com/anowacki/SeisSplit.jl) performs shear
  wave splitting analysis.
- [Beamforming](https://github.com/anowacki/Beamforming.jl) uses arrays
  of seismic recordings to stack data to detect signals.
- [SeisTau](https://github.com/anowacki/SeisTau.jl) calculates travel times
  of seismic waves through models of the Earth (and other planets).

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