# Job Training and Earnings

```@meta
CurrentModule = CausalGraphs
```

This vignette summarizes the National Supported Work (NSW) job-training
example. The full Quarto source is in `vignettes/Job_Training_NSW.qmd`.

## Causal Question

What is the effect of job-training participation on real earnings in 1978?

The treatment is `treat`, the outcome is `re78`, and the baseline covariates are
age, education, race/ethnicity indicators, marital status, degree status, and
pre-treatment real earnings in 1974 and 1975.

## Experimental Graph

If the NSW sample is treated as a randomized experiment, the graph is:

```text
A -> Y
```

```@example nsw_graph
using CausalGraphs

experimental_graph = make_graph(
    vertices = [:treat, :re78],
    di_edges = [(:treat, :re78)],
)

identify(experimental_graph, :treat, :re78).strategy
```

```@example nsw_graph
markov_pillow(experimental_graph, :treat; treatment=:treat)
```

No baseline adjustment is required under this graph.

## Observational Selection Graph

For an observational analysis, the graph allows measured baseline variables to
affect both program participation and later earnings:

```text
X -> A -> Y
X -----> Y
```

```@example nsw_graph
vertices = [
    :age, :educ, :black, :hisp, :marr, :nodegree,
    :re74, :re75, :treat, :re78,
]

baseline = setdiff(vertices, [:treat, :re78])
edges = vcat(
    [(w, :treat) for w in baseline],
    [(w, :re78) for w in baseline],
    [(:treat, :re78)],
)

observational_graph = make_graph(vertices=vertices, di_edges=edges)
identify(observational_graph, :treat, :re78).strategy
```

```@example nsw_graph
markov_pillow(observational_graph, :treat; treatment=:treat)
```

Under this graph, the effect is identified by adjusting for the measured
baseline covariates.

## Unmeasured Selection

If unmeasured motivation, caseworker discretion, health, networks, or local
labor-market opportunity confound training and later earnings, encode that as a
bidirected edge:

```text
X -> A -> Y
X -----> Y
A <-> Y
```

```@example nsw_graph
selection_graph = make_graph(
    vertices = vertices,
    di_edges = edges,
    bi_edges = [(:treat, :re78)],
)

identify(selection_graph, :treat, :re78).strategy
```

With that additional assumption, `CausalGraphs.jl` reports the effect as not
identified by the implemented criteria.

## Estimation

The full Quarto vignette downloads the NSW data from RDatasets and estimates the
effect under both identified graphs. In a local run, the experimental graph
gives a TMLE ACE of about 1794, while the baseline-adjusted observational graph
gives a TMLE ACE of about 1638. The key point is not the small numerical
difference in this sample, but that each estimate is tied to a different graph.

## Takeaway

The NSW example is useful because it makes the identifying assumptions visible:

- randomized assignment implies no adjustment;
- measured selection implies adjustment for baseline covariates;
- unmeasured selection makes the effect non-identified in this ADMG.
