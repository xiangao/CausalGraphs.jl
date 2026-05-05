# Nested-Fixable Estimation and Missing Data

```@meta
CurrentModule = CausalGraphs
```

This vignette shows the two pieces that go beyond standard backdoor/front-door
workflows: nested-fixable effects and missing-data weighting. All datasets are
simulated.

The common theme is identification before estimation. The graph determines
whether the target law can be represented using observed data. Only after that
step does the package construct weights and nuisance regressions.

## Nested-Fixable Effects

A treatment is nested-fixable when the target effect is identified by a
One-Line ID style fixing sequence, but the treatment is neither a-fixable nor
p-fixable. This happens in ADMGs where hidden-confounding structure is too
complex for ordinary backdoor adjustment or a front-door/NPS argument.

The word "nested" comes from nested Markov models for ADMGs. Identification is
described by repeatedly "fixing" variables: moving a variable from random to
fixed status and updating the kernel for the remaining random variables. If a
valid fixing sequence isolates the target counterfactual distribution, the
effect is identified.

```@example nested_missing
using CausalGraphs, DataFrames, Random

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
```

```@example nested_missing
draw_graph(graph; direction="TB")
```

For nested-fixable effects, `estimate_causal()` returns augmented nested IPW
(ANIPW) and nested IPW (NIPW) estimates.

Nested rebalancing weights are analogous to ordinary IPW, but they balance the
pieces of the observed law that appear in the nested identification formula.
ANIPW augments those weights with outcome regression terms, in the same spirit
that AIPW augments ordinary IPW.

```@example nested_missing
Random.seed!(42)
n = 500
U_EI = randn(n); U_ET = randn(n); U_ITox = randn(n)
U_VA = randn(n); U_VCD4 = randn(n)
ViralLoad = randn(n)
Exercise  = randn(n) .+ 0.5 .* U_EI .+ 0.5 .* U_ET
Income    = randn(n) .+ 0.4 .* ViralLoad .+ 0.4 .* U_EI .+ 0.4 .* U_ITox
T         = Float64.(rand(n) .< 1 ./(1 .+ exp.(-(0.5 .* ViralLoad .+ 0.3 .* Income .+ 0.4 .* U_ET))))
Toxicity  = randn(n) .+ 0.6 .* T .+ 0.3 .* U_ITox
Adherence = randn(n) .+ 0.5 .* Toxicity .+ 0.3 .* U_VA
CD4       = randn(n) .+ 0.4 .* T .+ 0.3 .* Adherence .+
            0.2 .* Exercise .+ 0.3 .* ViralLoad .+ 0.3 .* U_VCD4 .+ 0.2 .* U_VA
df = DataFrame(ViralLoad=ViralLoad, Income=Income, Exercise=Exercise,
               T=T, Toxicity=Toxicity, Adherence=Adherence, CD4=CD4)

res = estimate_causal(a=[1, 0], data=df, graph=graph, treatment=:T, outcome=:CD4)
r = x -> round(x, sigdigits=4)
(ACE=r(res[:ANIPW].ACE), lower_ci=r(res[:ANIPW].lower_ci), upper_ci=r(res[:ANIPW].upper_ci))
```

## Missing-Data Weighting

An mDAG extends a causal graph with missingness indicators. A variable such as
`X` has an indicator `Rx` showing whether `X` is observed. Edges into `Rx`
encode the missingness mechanism.

The key identification question is whether the target full-data law can be
recovered from the observed-data law. If it can, the package estimates
observation probabilities and forms inverse-probability weights. Those weights
can be passed to `estimate_causal(...; sample_weights=wts)` so complete cases
are weighted back toward the target population.

This separates two assumptions: the causal ADMG describes treatment, outcome,
and covariates; the mDAG describes how values become missing.

```@example nested_missing
mdag = make_mdag(
    obs_variables = ["A", "Y"],
    missing_variables = ["X"],
    missing_indicators = ["Rx"],
    di_edges = [("X", "A"), ("X", "Y"), ("X", "Rx"), ("A", "Y")],
)

ID = ID_algorithm(mdag)
ID.target.id_status
```

```@example nested_missing
draw_graph(mdag)
```

`compute_missing_weights()` estimates the identified missingness propensities
and returns inverse-probability weights that can be passed to
`estimate_causal(...; sample_weights=wts)`.
