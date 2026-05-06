# CausalGraphs.jl

`CausalGraphs.jl` is a Julia package for graph-based causal identification,
semiparametric effect estimation, and missing-data weighting in acyclic directed
mixed graphs (ADMGs) and missingness DAGs (mDAGs).
It includes estimator-routing checks for a-fixable, p-fixable, and
nested-fixable effects, plus symbolic Pearl-Shpitser ID for more general ADMG
queries.

It brings together ideas and workflows from Anna Guo and Razieh Nabi's R
packages [`flexCausal`](https://github.com/annaguo-bios/flexCausal) and
[`flexMissing`](https://github.com/annaguo-bios/flexMissing), rewritten and
adapted for Julia.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/xiangao/CausalGraphs.jl")
```

## Quick Start

```@example causalgraphs_home
using CausalGraphs, DataFrames, Random

Random.seed!(1)
n = 500
X = randn(n)
A = Float64.(rand(n) .< 1 ./ (1 .+ exp.(-X)))
Y = 2 .* A .+ X .+ randn(n)
data = DataFrame(X=X, A=A, Y=Y)

graph = make_graph(
    vertices = [:X, :A, :Y],
    di_edges = [(:X, :A), (:X, :Y), (:A, :Y)],
)

identify(graph, :A, :Y).strategy
```

```@example causalgraphs_home
result = estimate_causal(
    a = [1, 0],
    data = data,
    graph = graph,
    treatment = :A,
    outcome = :Y,
)

result[:TMLE].ACE
```

## Vignettes

| Vignette | Description |
|----------|-------------|
| [Graphs and Identification](vignettes/Graphs_and_Identification.md) | ADMG construction, visualization, graph properties, and `identify()` routing |
| [General ADMG ID Algorithm](vignettes/ADMG_ID_Algorithm.md) | Symbolic Pearl-Shpitser ID, hedge failures, fixing sequences, and nested-fixability diagnostics |
| [Estimation: Backdoor and Front-Door](vignettes/Estimation_Backdoor_Frontdoor.md) | Backdoor/a-fixable and front-door/p-fixable estimation workflows |
| [Nested-Fixable and Missing Data](vignettes/Nested_and_Missing_Data.md) | Nested-fixable effects, ANIPW/NIPW, and missing-data weighting with mDAGs |
| [Smoking Cessation NHEFS](vignettes/Smoking_Cessation_NHEFS.md) | End-to-end real-data workflow: hypothesize a graph, identify, estimate, and compare assumptions |
| [Job Training NSW](vignettes/Job_Training_NSW.md) | Economics example: experimental assignment, measured selection, unmeasured selection, identification, and estimation |
| [Job Search Mediation JOBS II](vignettes/Job_Search_Mediation_JOBS.md) | Economics mediation example: DAG, total-effect estimation, and natural direct and indirect effects |

## Core API

- `make_graph`
- `draw_graph`
- `to_mermaid`
- `identify`
- `ID_algorithm`
- `estimate_causal`
- `make_mdag`
- `compute_missing_weights`

## References

- Anna Guo and Razieh Nabi,
  [`flexCausal`](https://github.com/annaguo-bios/flexCausal): causal effect
  estimation in ADMGs with hidden variables.
- Anna Guo and Razieh Nabi,
  [`flexMissing`](https://github.com/annaguo-bios/flexMissing):
  weighting-based identification and estimation in graphical models of missing
  data.
