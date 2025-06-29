using TimeseriesUtilities
using Documenter

DocMeta.setdocmeta!(TimeseriesUtilities, :DocTestSetup, :(using TimeseriesUtilities); recursive=true)

makedocs(;
    modules=[TimeseriesUtilities],
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
            "outliers.md"
        ]
    ],
)

deploydocs(;
    repo="github.com/Beforerr/TimeseriesUtilities.jl",
    push_preview = true
)
