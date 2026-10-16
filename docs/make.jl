using TimeseriesUtilities
using Documenter

DocMeta.setdocmeta!(TimeseriesUtilities, :DocTestSetup, :(using TimeseriesUtilities); recursive=true)

# Home is README.md with its walkthrough blocks evaluated, so the README stays the one source.
let readme = read(joinpath(pkgdir(TimeseriesUtilities), "README.md"), String)
    parts = split(readme, "\n## Walkthrough\n"; limit = 2)
    length(parts) == 2 || error("README.md has no `## Walkthrough` section to evaluate")
    walkthrough = replace(parts[2], "```julia\n" => "```@repl walkthrough\n")
    write(
        joinpath(@__DIR__, "src", "index.md"),
        "```@meta\nEditURL = \"../../README.md\"\n```\n\n", parts[1], "\n## Walkthrough\n", walkthrough,
    )
end

# Documenter's `@repl` output is truncated to `displaysize()`, which reads LINES; keep long arrays short.
ENV["LINES"] = "12"

makedocs(;
    modules=[TimeseriesUtilities],
    # `@autodocs` renders every docstring but the module one, which is README.md, rendered as Home
    checkdocs=:none,
    authors="Beforerr <zzj956959688@gmail.com> and contributors",
    sitename="TimeseriesUtilities.jl",
    format=Documenter.HTML(;
        canonical="https://Beforerr.github.io/TimeseriesUtilities.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => [
            "outliers.md",
            "cadence.md"
        ],
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/Beforerr/TimeseriesUtilities.jl",
    push_preview = true
)
