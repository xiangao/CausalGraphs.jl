# Nested-Fixable Estimation and Missing Data

```@meta
CurrentModule = CausalGraphs
```

This vignette shows the two pieces that go beyond standard backdoor/front-door
workflows: nested-fixable effects and missing-data weighting.

## Nested-Fixable Effects

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

For nested-fixable effects, `estimate_causal()` returns nested IPW and augmented
nested IPW estimates.

## Missing-Data Weighting

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
