using Documenter, Seis

makedocs(
    sitename = "Seis.jl documentation",
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md",
        "Plotting with Plots.jl" => "plotting-plots.md",
        "Plotting with Makie.jl" => "plotting-makie.md",
        "Function index" => "function-index.md",
        "Internals" => "internal-index.md",
    ],
    # Workaround problem including `add_picks!` docstring
    # TODO: Remove this if and when SeisTau becomes a module extension
    remotes = nothing,
    format = Documenter.HTML(repolink="https://github.com/anowacki/Seis.jl.git"),
)

deploydocs(
    repo = "github.com/anowacki/Seis.jl.git",
)
