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

```@example nhefs_graph
draw_graph(conceptual_graph)
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

```@example nhefs_graph
draw_graph(graph; direction="TB")
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

```@example nhefs_graph
draw_graph(sensitivity_graph; direction="TB")
```

The purpose of the DAG is to make this assumption explicit before estimation.

## Estimation

```@example nhefs_graph
using DataFrames, Random

Random.seed!(1)
n = 1566
sex   = rand(0:1, n)
race  = rand(0:1, n)
age   = rand(25:74, n)
school        = rand(6:17, n)
smokeintensity = rand(1:40, n)
smokeyrs      = rand(1:40, n)
exercise      = rand(0:2, n)
active        = rand(0:2, n)
wt71          = 55 .+ 20 .* randn(n)
qsmk = Float64.(rand(n) .< 1 ./(1 .+ exp.(-(
    -1.5 .+ 0.01 .* age .- 0.01 .* smokeintensity .- 0.01 .* smokeyrs))))
wt82_71 = 2.5 .* qsmk .+ 0.02 .* (age .- 45) .+ 0.03 .* smokeintensity .+ randn(n) .* 5

data = DataFrame(sex=sex, race=race, age=age, school=school,
                 smokeintensity=smokeintensity, smokeyrs=smokeyrs,
                 exercise=exercise, active=active, wt71=wt71,
                 qsmk=qsmk, wt82_71=wt82_71)

res = estimate_causal(a=[1, 0], data=data, graph=graph,
                      treatment=:qsmk, outcome=:wt82_71)
(ACE=res[:TMLE].ACE, lower_ci=res[:TMLE].lower_ci, upper_ci=res[:TMLE].upper_ci)
```

The true ACE in this simulation is 2.5 kg. On the real NHEFS data the
TMLE estimate is approximately 3.5 kg (95% CI 2.4–4.5).

