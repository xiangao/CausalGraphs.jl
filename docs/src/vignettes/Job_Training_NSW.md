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
draw_graph(experimental_graph; direction="LR")
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
conceptual_graph = make_graph(
    vertices = [:X, :A, :Y],
    di_edges = [(:X, :A), (:X, :Y), (:A, :Y)],
)

draw_graph(conceptual_graph; direction="LR")
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
draw_graph(observational_graph; direction="LR")
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

```@example nsw_graph
draw_graph(selection_graph; direction="LR")
```

With that additional assumption, `CausalGraphs.jl` reports the effect as not
identified by the implemented criteria.

## Estimation

The NSW data can be loaded directly from RDatasets. The estimate under the
experimental graph uses only `treat` and `re78`, while the observational graph
also adjusts for the baseline covariates through the Markov pillow of `treat`.

```@example nsw_graph
using DataFrames, DelimitedFiles, Downloads, Statistics

function read_rdatasets_csv(url, cols)
    local_file = Downloads.download(url; timeout=120)
    x, header = readdlm(local_file, ',', Any, '\n'; header=true)
    raw = DataFrame(x, Symbol.(vec(header)))
    DataFrame((c => Float64.(raw[!, c]) for c in cols)...)
end

url = "https://vincentarelbundock.github.io/Rdatasets/csv/causaldata/nsw_mixtape.csv"
cols = [
    :treat, :age, :educ, :black, :hisp, :marr, :nodegree,
    :re74, :re75, :re78,
]

data = read_rdatasets_csv(url, cols)

(n = nrow(data),
 treated = Int(sum(data.treat .== 1)),
 controls = Int(sum(data.treat .== 0)))
```

Under the experimental graph, no baseline adjustment is made.

```@example nsw_graph
raw_difference = mean(data.re78[data.treat .== 1]) -
                 mean(data.re78[data.treat .== 0])

experimental_res = estimate_causal(
    a = [1, 0],
    data = select(data, [:treat, :re78]),
    graph = experimental_graph,
    treatment = :treat,
    outcome = :re78,
)

r = x -> round(x, sigdigits=4)

(raw_difference = r(raw_difference),
 TMLE_ACE = r(experimental_res[:TMLE].ACE),
 lower_ci = r(experimental_res[:TMLE].lower_ci),
 upper_ci = r(experimental_res[:TMLE].upper_ci))
```

Under the observational selection-on-observables graph, the estimator adjusts
for the measured baseline variables.

```@example nsw_graph
observational_res = estimate_causal(
    a = [1, 0],
    data = data,
    graph = observational_graph,
    treatment = :treat,
    outcome = :re78,
)

(TMLE_ACE = r(observational_res[:TMLE].ACE),
 lower_ci = r(observational_res[:TMLE].lower_ci),
 upper_ci = r(observational_res[:TMLE].upper_ci),
 Onestep_ACE = r(observational_res[:Onestep].ACE),
 Gcomp_ACE = r(observational_res[:Gcomp].ACE),
 IPW_ACE = r(observational_res[:IPW].ACE))
```

The key point is not the small numerical difference in this sample. It is that
each estimate is tied to a different graph. Under the unmeasured-selection graph
above, the effect is not identified, so `estimate_causal` should not be run.

## Takeaway

The NSW example is useful because it makes the identifying assumptions visible:

- randomized assignment implies no adjustment;
- measured selection implies adjustment for baseline covariates;
- unmeasured selection makes the effect non-identified in this ADMG.
