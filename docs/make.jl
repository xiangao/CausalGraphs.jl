using Documenter
using CausalGraphs

makedocs(
    sitename = "CausalGraphs.jl",
    modules = [CausalGraphs],
    format = Documenter.HTML(
        assets = [Documenter.asset(
            "https://cdn.jsdelivr.net/npm/mermaid@11/dist/mermaid.min.js",
            class = :js, islocal = false,
        )],
    ),
    pages = [
        "Home" => "index.md",
        "Vignettes" => [
            "Graphs and Identification" => "vignettes/Graphs_and_Identification.md",
            "Estimation: Backdoor and Front-Door" => "vignettes/Estimation_Backdoor_Frontdoor.md",
            "Nested-Fixable and Missing Data" => "vignettes/Nested_and_Missing_Data.md",
            "Smoking Cessation NHEFS" => "vignettes/Smoking_Cessation_NHEFS.md",
        ],
        "Reference" => "reference.md",
    ],
    warnonly = true,
)

deploydocs(
    repo = "github.com/xiangao/CausalGraphs.jl.git",
    devbranch = "main",
    push_preview = false,
)
