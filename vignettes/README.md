# CausalGraphs.jl Vignettes

Full documentation: **https://xiangao.github.io/CausalGraphs.jl/dev/**

| Vignette | Description |
|---|---|
| [Graphs and Identification](Graphs_and_Identification.html) | ADMG construction, graph visualization, graph properties, and `identify()` routing |
| [General ADMG ID Algorithm](ADMG_ID_Algorithm.html) | Symbolic Pearl-Shpitser ID, finite-support plug-in estimation, hedge failures, and fixing diagnostics |
| [Real Example: Berkeley Admissions with ID](Berkeley_Admissions_ID.html) | Real-data front-door ADMG, Pearl-Shpitser ID functional, finite-support plug-in estimation, and non-identification sensitivity |
| [Estimation: Backdoor and Front-Door](Estimation_Backdoor_Frontdoor.html) | Backdoor/a-fixable and front-door/p-fixable estimation workflows |
| [Nested-Fixable Estimation and Missing Data](Nested_and_Missing_Data.html) | Nested-fixable effects, ANIPW/NIPW, and missing-data weighting with mDAGs |
| [Real Example: Smoking Cessation and Weight Change](Smoking_Cessation_NHEFS.html) | End-to-end NHEFS example: hypothesize a graph, identify, estimate, and compare assumptions |
| [Economics Example: Job Training and Earnings](Job_Training_NSW.html) | NSW job-training example: experimental assignment, measured selection, unmeasured selection, identification, and estimation |
| [Economics Example: Job Search Mediation](Job_Search_Mediation_JOBS.html) | JOBS II mediation example: DAG, total-effect identification and estimation, natural direct and indirect effects |

The `.qmd` files in this directory are the Quarto sources for the rendered HTML vignettes.
