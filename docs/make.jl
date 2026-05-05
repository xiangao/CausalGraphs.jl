using Documenter
using CausalGraphs

makedocs(
    sitename = "CausalGraphs.jl",
    modules = [CausalGraphs],
    pages = [
        "Home" => "index.md",
        "Vignettes" => [
            "Graphs and Identification" => "vignettes/Graphs_and_Identification.md",
            "Estimation: Backdoor and Front-Door" => "vignettes/Estimation_Backdoor_Frontdoor.md",
            "Nested-Fixable and Missing Data" => "vignettes/Nested_and_Missing_Data.md",
            "Smoking Cessation NHEFS" => "vignettes/Smoking_Cessation_NHEFS.md",
            "Job Training NSW" => "vignettes/Job_Training_NSW.md",
            "Job Search Mediation JOBS II" => "vignettes/Job_Search_Mediation_JOBS.md",
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
