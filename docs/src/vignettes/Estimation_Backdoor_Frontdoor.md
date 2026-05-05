# Estimation: Backdoor and Front-Door

```@meta
CurrentModule = CausalGraphs
```

This vignette covers the common `estimate_causal()` workflows for backdoor and
front-door style graphs. All datasets are simulated with known true effects.

## Backdoor

```@example backdoor_frontdoor
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

```@example backdoor_frontdoor
draw_graph(graph)
```

```@example backdoor_frontdoor
res = estimate_causal(
    a = [1, 0],
    data = data,
    graph = graph,
    treatment = :A,
    outcome = :Y,
)

r = x -> round(x, sigdigits=4)
(ACE=r(res[:TMLE].ACE), lower_ci=r(res[:TMLE].lower_ci), upper_ci=r(res[:TMLE].upper_ci))
```

The true ACE is 2.0. For an a-fixable treatment, `estimate_causal()` uses the
Markov pillow of the treatment as the adjustment set and returns TMLE,
one-step/AIPW, G-computation, and IPW estimates.

## Front-Door / P-Fixable

```@example backdoor_frontdoor
Random.seed!(2)
n = 700
A = Float64.(rand(n) .< 0.5)
M = 0.8 .* A .+ randn(n)
Y = 1.5 .* M .+ randn(n)
fd_data = DataFrame(A=A, M=M, Y=Y)

fd_graph = make_graph(
    vertices = [:A, :M, :Y],
    di_edges = [(:A, :M), (:M, :Y)],
    bi_edges = [(:A, :Y)],
)

identify(fd_graph, :A, :Y).strategy
```

```@example backdoor_frontdoor
draw_graph(fd_graph)
```

```@example backdoor_frontdoor
fd_res = estimate_causal(
    a = [1, 0],
    data = fd_data,
    graph = fd_graph,
    treatment = :A,
    outcome = :Y,
)

r = x -> round(x, sigdigits=4)
(ACE=r(fd_res[:TMLE].ACE), lower_ci=r(fd_res[:TMLE].lower_ci), upper_ci=r(fd_res[:TMLE].upper_ci))
```

The true ACE is 1.2 (= 0.8 × 1.5, by the front-door formula). The package
routes to the NPS/front-door style TMLE for p-fixable graphs.

