"""
# `SAC`

Module for dealing with SAC data.
"""
module SAC

using Statistics: mean

import Glob

import ..Seis
import Geodesics

include("constants.jl")
include("types.jl")
include("io.jl")
include("operations.jl")

end # module
