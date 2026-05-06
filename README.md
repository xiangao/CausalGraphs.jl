# CausalGraphs.jl

`CausalGraphs.jl` is a Julia package for identification and estimation of
causal effects in acyclic directed mixed graphs (ADMGs), with support for hidden
variables, nested-fixable effects, and missing-data weighting.

It combines graph-based causal identification, semiparametric effect
estimation, and missing-data weighting in one package:

1. Build an ADMG with directed and bidirected edges.
2. Check whether an effect is a-fixable, p-fixable, nested-fixable, identified
   by the general ID algorithm, or not identified.
3. Route automatically to the appropriate semiparametric estimator.
4. Optionally identify missing-data mechanisms on an mDAG and pass
   inverse-probability weights into the causal estimator.

The package currently implements graph utilities, automatic identification
routing, backdoor TMLE, front-door/NPS TMLE, nested ANIPW/NIPW, default
parametric nuisance fits, MLJ-based SuperLearner ensembles, symbolic
Pearl-Shpitser ID for ADMGs, finite-support plug-in estimation for discrete ID
functionals, mDAG missing-data identification, missingness propensity
estimation, and missing-data weighting.

`CausalGraphs.jl` brings together ideas and workflows from Anna Guo and Razieh
Nabi's R packages [`flexCausal`](https://github.com/annaguo-bios/flexCausal)
and [`flexMissing`](https://github.com/annaguo-bios/flexMissing), rewritten and
adapted for Julia.

## Install/Load Locally

From this directory:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using CausalGraphs
```

Run tests:

```julia
using Pkg
Pkg.test()
```

## Vignettes

Full documentation: **https://xiangao.github.io/CausalGraphs.jl/dev/**

| Vignette | Description |
|---|---|
| [Graphs and Identification](https://xiangao.github.io/CausalGraphs.jl/dev/vignettes/Graphs_and_Identification/) | ADMG construction, graph visualization, graph properties, and `identify()` routing |
| [General ADMG ID Algorithm](https://xiangao.github.io/CausalGraphs.jl/dev/vignettes/ADMG_ID_Algorithm/) | Symbolic Pearl-Shpitser ID, finite-support plug-in estimation, hedge failures, and fixing diagnostics |
| [Estimation: Backdoor and Front-Door](https://xiangao.github.io/CausalGraphs.jl/dev/vignettes/Estimation_Backdoor_Frontdoor/) | Backdoor/a-fixable and front-door/p-fixable estimation workflows |
| [Nested-Fixable Estimation and Missing Data](https://xiangao.github.io/CausalGraphs.jl/dev/vignettes/Nested_and_Missing_Data/) | Nested-fixable effects, ANIPW/NIPW, and missing-data weighting with mDAGs |
| [Real Example: Smoking Cessation and Weight Change](https://xiangao.github.io/CausalGraphs.jl/dev/vignettes/Smoking_Cessation_NHEFS/) | End-to-end NHEFS example: hypothesize a graph, identify, estimate, and compare assumptions |
| [Economics Example: Job Training and Earnings](https://xiangao.github.io/CausalGraphs.jl/dev/vignettes/Job_Training_NSW/) | NSW job-training example: experimental assignment, measured selection, unmeasured selection, identification, and estimation |
| [Economics Example: Job Search Mediation](https://xiangao.github.io/CausalGraphs.jl/dev/vignettes/Job_Search_Mediation_JOBS/) | JOBS II mediation example: DAG, total-effect identification and estimation, natural direct and indirect effects |

The Quarto sources live in [`vignettes/`](vignettes/README.md).

## Quick Start

```julia
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

draw_graph(graph)

id = identify(graph, :A, :Y)
id.strategy

result = estimate_causal(
    a = [1, 0],
    data = data,
    graph = graph,
    treatment = :A,
    outcome = :Y,
)

result[:TMLE].ACE
result[:TMLE].lower_ci
result[:TMLE].upper_ci
```

## Features

### Graphs

- `make_graph()` constructor for ADMGs
- Mermaid and Graphviz DOT graph export with `draw_graph()`, `to_mermaid()`,
  and `to_dot()`
- Directed and bidirected edges
- Ancestors, descendants, parents, children
- Districts and all districts
- Markov blankets and Markov pillows
- Subgraphs and fixed graphs
- Reachable closures
- Fixability and primal-fixability checks
- General fixing sequences
- Nonparametric saturation and Markov-blanket shielding checks

### Identification

`identify(graph, treatment, outcome)` returns a named tuple whose `strategy` is
one of:

| Strategy | Meaning |
|---|---|
| `:a_fixable` | Backdoor/a-fixable effect |
| `:p_fixable` | Front-door, primal-fixable, or NPS effect |
| `:nested_fixable` | Identified by a One-Line-ID style nested check |
| `:id_algorithm` | Identified by the general Pearl-Shpitser ID algorithm; routed to finite-support plug-in estimation for discrete variables |
| `:not_identified` | Not identified by the implemented criteria |

For symbolic identification without estimator routing, use
`ID_algorithm(graph, treatment, outcome)`. For direct finite-support plug-in
estimation of a symbolic ID functional, use `estimate_id(...)`.

### Estimation

The main entry point is `estimate_causal(...)`. It routes automatically based on
the graph:

| Strategy | Returned estimators |
|---|---|
| a-fixable | `:TMLE`, `:Onestep`, `:Gcomp`, `:IPW` |
| p-fixable | `:TMLE`, `:Onestep` |
| nested-fixable | `:ANIPW`, `:NIPW` |
| ID-algorithm | `:IDPlugin` for finite-support/discrete variables |

The ID plug-in estimator enumerates observed supports, so it is appropriate for
discrete variables. It is not a general continuous-variable TMLE.

The argument `a` can be a scalar, such as `a=1`, for `E[Y(a)]`, or a length-two
vector, such as `a=[1,0]`, for an ACE contrast.

Treatments are currently expected to be binary and coded as `0/1`.

### Nuisance Models

Default nuisance fits are implemented directly in Julia for numeric
`DataFrame` columns:

- linear regression for continuous outcomes,
- logistic regression for binary outcomes and propensity scores,
- weighted GLM/ridge-style nuisance fitting when `sample_weights` are supplied.

MLJ-based stacked ensembles are available through:

```julia
superlearner(binary=true)
superlearner(binary=false)
```

Cross-fitting is available for MLJ/SuperLearner fits. `sample_weights` are
currently supported for the default/GLM nuisance fits, but not for cross-fitted
MLJ models.

### Missing Data

`make_mdag()`, `ID_algorithm()`, and `propensity()` implement missing-data
identification and missingness propensity estimation inside CausalGraphs.jl.
`compute_missing_weights()` converts those propensities into inverse-probability
weights that can be passed to `estimate_causal`.

```julia
using CausalGraphs, DataFrames

mdag = make_mdag(
    obs_variables = ["A", "Y"],
    missing_variables = ["X"],
    missing_indicators = ["Rx"],
    di_edges = [("X", "A"), ("X", "Y"), ("X", "Rx"), ("A", "Y")],
)

ID = ID_algorithm(mdag)
wts = compute_missing_weights(mdag, data_with_missing;
                              ID=ID, complete_cases_only=true)

result = estimate_causal(
    a = [1, 0],
    data = dropmissing(data_with_missing),
    graph = graph,
    treatment = :A,
    outcome = :Y,
    sample_weights = wts,
)
```

## Examples

### Front-Door / NPS

```julia
graph = make_graph(
    vertices = [:X, :A, :M, :Y],
    bi_edges = [(:A, :Y)],
    di_edges = [(:X, :A), (:X, :Y), (:X, :M), (:A, :M), (:M, :Y)],
)

identify(graph, :A, :Y).strategy

result = estimate_causal(
    a = [1, 0],
    data = data,
    graph = graph,
    treatment = :A,
    outcome = :Y,
)

result[:TMLE].ACE
```

### Nested-Fixable Effect

```julia
graph = make_graph(
    vertices = [:ViralLoad, :Income, :Exercise, :T, :Toxicity, :Adherence, :CD4],
    di_edges = [
        (:ViralLoad, :T), (:ViralLoad, :CD4),
        (:Income, :T), (:Exercise, :CD4),
        (:T, :Toxicity), (:T, :CD4),
        (:Toxicity, :Adherence), (:Adherence, :CD4),
    ],
    bi_edges = [
        (:Income, :Toxicity),
        (:Exercise, :Income),
        (:Exercise, :T),
        (:ViralLoad, :Adherence),
        (:ViralLoad, :CD4),
    ],
)

identify(graph, :T, :CD4).strategy

result = estimate_causal(
    a = [1, 0],
    data = data,
    graph = graph,
    treatment = :T,
    outcome = :CD4,
)

result[:ANIPW].ACE
```

## How The Method Works

The graph determines the identification strategy before estimation starts.

For a-fixable effects, the package uses the Markov pillow of the treatment as
the adjustment set. It estimates `E[Y | A, mp(A)]` and `P(A | mp(A))`, then
returns plug-in/G-computation, IPW, one-step/AIPW, and TMLE estimates.

For p-fixable effects, the graph is decomposed into pre-treatment variables,
variables in the treatment district, and the remaining post-treatment
variables. The NPS estimator fits an outcome regression and sequential
regressions along a topological order, using density-ratio terms to construct
the observed-data estimating equation. The TMLE iteratively targets the
propensity, outcome regression, and sequential regressions until the empirical
estimating equation is small or the iteration limit is reached.

For nested-fixable effects, the package uses a One-Line-ID style check to
construct the relevant outcome ancestors and modified districts. It computes
nested rebalancing weights through intrinsic kernels, then estimates the effect
with nested IPW and augmented nested IPW.

Missing-data adjustment is handled by the mDAG tools in CausalGraphs.jl:
the package identifies missingness propensities, estimates them, converts them
to row weights, and passes those weights through the causal estimator.

## R-to-Julia Naming

The package follows Julia-style snake_case names:

| Concept | Julia |
|---|---|
| ADMG constructor | `make_graph` |
| One-Line ID style routing | `identify` |
| Causal effect estimation | `estimate_causal` |
| Markov blanket | `markov_blanket` |
| Markov pillow | `markov_pillow` |
| mDAG constructor | `make_mdag` |
| Missing-data identification | `ID_algorithm` |
| Missing-data IPW | `compute_missing_weights` |
| SuperLearner stack | `superlearner` |

## Current Limitations

- Treatments are currently binary and coded as `0/1`.
- Kernel `densratio` is not implemented.
- `sample_weights` are supported for default/GLM nuisance fits, but not for
  cross-fitted MLJ/SuperLearner nuisance fits.

## References

The implementation is based on ideas from:

- Bhattacharya, Nabi, and Shpitser. *Semiparametric Inference for Causal Effects
  in Graphical Models with Hidden Variables*.
- Shpitser and Pearl. *Identification of Joint Interventional Distributions in
  Recursive Semi-Markovian Causal Models*.
- Richardson, Robins, and Shpitser. *Nested Markov Properties for Acyclic
  Directed Mixed Graphs*.
- Guo and Nabi. *Average Causal Effect Estimation in DAGs with Hidden
  Variables: Extensions of Back-Door and Front-Door Criteria*.
- Guo and Nabi. *Weighting-Based Identification and Estimation in Graphical
  Models of Missing Data*.
- Anna Guo and Razieh Nabi's R package
  [`flexCausal`](https://github.com/annaguo-bios/flexCausal), for causal effect
  estimation in ADMGs with hidden variables.
- Anna Guo and Razieh Nabi's R package
  [`flexMissing`](https://github.com/annaguo-bios/flexMissing), for
  weighting-based identification and estimation in graphical models of missing
  data.
