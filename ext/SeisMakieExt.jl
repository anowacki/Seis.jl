"""
# `SeisMakieExt`

Extension module providing plots using Makie when a Makie backend
module is loaded.
"""
module SeisMakieExt

import Dates
import Seis

@static if isdefined(Base, :get_extension)
    import Makie
else
    import ..Makie
end

include("SeisMakieExt/util.jl")
include("SeisMakieExt/plot_traces.jl")
include("SeisMakieExt/plot_section.jl")
include("SeisMakieExt/plot_hodogram.jl")
include("SeisMakieExt/pick_plot.jl")

end # module
