# Estimation: Backdoor and Front-Door

```@meta
CurrentModule = CausalGraphs
```

This vignette covers the common `estimate_causal()` workflows for backdoor and
front-door style graphs. All datasets are simulated with known true effects.

The workflow is: encode assumptions in an ADMG, call `identify()`, then call
`estimate_causal()` only when the graph identifies the target effect. The
estimand here is the average causal effect `E[Y(1)-Y(0)]`.

## Backdoor

Backdoor identification applies when observed pre-treatment variables block all
noncausal paths from treatment to outcome. In this graph, `X` affects both `A`
and `Y`, so a raw treated-control comparison is confounded. The adjustment
formula is:

```text
E[Y(a)] = sum_x E[Y | A=a, X=x] P(X=x)
```

The a-fixable check is the ADMG version of this adjustment logic. The package
uses the treatment's Markov pillow as the adjustment set.

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

G-computation models the outcome regression. IPW models the propensity score.
One-step/AIPW and TMLE combine both nuisance models; TMLE additionally updates
the outcome regression to solve the efficient influence-function equation.

## Front-Door / P-Fixable

Front-door identification is useful when treatment and outcome are confounded,
but the causal effect flows through observed post-treatment variables. The
bidirected edge `A <-> Y` means ordinary backdoor adjustment fails. The
mediator `M` helps because the graph identifies how `A` changes `M` and how
`M` changes `Y`.

The classical front-door functional is:

```text
E[Y(a)] = sum_m P(M=m | A=a) sum_a' E[Y | M=m, A=a'] P(A=a')
```

The package uses the broader p-fixable/NPS route. A treatment is p-fixable when
none of its children are in its bidirected-connected district. Intuitively,
treatment may be confounded with the outcome, but its immediate causal children
are not hidden-confounded with treatment.

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
