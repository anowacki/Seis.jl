using Documenter, Seis

makedocs(
    sitename = "Seis.jl documentation",
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md",
        "Plotting" => "plotting.md",
        "Function index" => "function-index.md",
        "Internals" => "internal-index.md",
        ]
    )

deploydocs(
    repo = "github.com/anowacki/Seis.jl.git",
)
