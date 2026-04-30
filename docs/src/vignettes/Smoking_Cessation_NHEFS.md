# Smoking Cessation and Weight Change

```@meta
CurrentModule = CausalGraphs
```

This vignette summarizes the NHEFS smoking-cessation example. The full Quarto
source is in `vignettes/Smoking_Cessation_NHEFS.qmd`.

## Causal Question

What is the effect of quitting smoking on weight change from 1971 to 1982?

The outcome is `wt82_71`, the treatment is `qsmk`, and the baseline covariates
include demographics, education, smoking history, activity, exercise, and
baseline weight.

## Conceptual Graph

The main adjustment graph is:

```text
X -> A -> Y
X -----> Y
```

where `X` is the baseline covariate set, `A` is quitting smoking, and `Y` is
later weight change. This graph is an assumption, not a fact learned from the
data.

```@example nhefs_graph
using CausalGraphs

conceptual_graph = make_graph(
    vertices = [:X, :A, :Y],
    di_edges = [(:X, :A), (:X, :Y), (:A, :Y)],
)

identify(conceptual_graph, :A, :Y).strategy
```

## Expanded Graph

For estimation, `X` is expanded into actual data columns:

```@example nhefs_graph
vertices = [
    :sex, :race, :age, :school,
    :smokeintensity, :smokeyrs,
    :exercise, :active, :wt71,
    :qsmk, :wt82_71,
]

baseline = setdiff(vertices, [:qsmk, :wt82_71])
edges = vcat(
    [(w, :qsmk) for w in baseline],
    [(w, :wt82_71) for w in baseline],
    [(:qsmk, :wt82_71)],
)

graph = make_graph(vertices=vertices, di_edges=edges)
identify(graph, :qsmk, :wt82_71).strategy
```

## Sensitivity Graph

If we add an unmeasured common cause of quitting and weight change, the effect is
not identified by the implemented criteria:

```@example nhefs_graph
sensitivity_graph = make_graph(
    vertices = vertices,
    di_edges = edges,
    bi_edges = [(:qsmk, :wt82_71)],
)

identify(sensitivity_graph, :qsmk, :wt82_71).strategy
```

The purpose of the DAG is to make this assumption explicit before estimation.

