# Graphs and Identification

```@meta
CurrentModule = CausalGraphs
```

This vignette introduces the graph layer and the `identify()` function.

A graph is a compact statement of causal assumptions. It tells us which
variables can directly affect other variables and which observed variables may
share hidden common causes. Identification asks whether an interventional
quantity, such as `E[Y(a)]` or `E[Y(1)-Y(0)]`, can be written in terms of the
observed data distribution under those assumptions.

## Constructing an ADMG

An acyclic directed mixed graph has directed edges for observed causal
relationships and bidirected edges for latent common causes.

In an ADMG:

- `A -> B` means `A` may directly cause `B`;
- `A <-> B` means `A` and `B` share an unmeasured common cause;
- directed arrows cannot form cycles.

ADMGs are useful because they represent latent confounding without drawing each
unobserved variable. Instead of adding a hidden `U` with `U -> A` and `U -> Y`,
we draw `A <-> Y`.

```@example graphs_identification
using CausalGraphs

g_fd = make_graph(
    vertices = [:A, :M, :Y],
    di_edges = [(:A, :M), (:M, :Y)],
    bi_edges = [(:A, :Y)],
)

draw_graph(g_fd; direction="LR")
```

```@example graphs_identification
g_bd = make_graph(
    vertices = [:X, :A, :Y],
    di_edges = [(:X, :A), (:X, :Y), (:A, :Y)],
)

draw_graph(g_bd; direction="LR")
```

## Graph Properties

Graph properties are the vocabulary used by identification algorithms. Parents
and children describe direct arrows. Ancestors and descendants describe all
directed upstream or downstream variables. Districts are bidirected-connected
components and encode latent confounding. Markov pillows are local adjustment
sets used by the a-fixable/backdoor estimator.

```@example graphs_identification
parents(g_fd, :M)
```

```@example graphs_identification
children(g_fd, :A)
```

```@example graphs_identification
district(g_fd, :A)
```

```@example graphs_identification
markov_pillow(g_bd, :A; treatment=:A)
```

## Identification

`identify()` routes the graph to the implemented estimation strategy.

```@example graphs_identification
identify(g_bd, :A, :Y).strategy
```

```@example graphs_identification
identify(g_fd, :A, :Y).strategy
```

The implemented strategies are:

| Strategy | Interpretation |
|----------|----------------|
| `:a_fixable` | Backdoor/a-fixable effect |
| `:p_fixable` | Front-door, primal-fixable, or NPS effect |
| `:nested_fixable` | Identified by a One-Line-ID style nested check |
| `:id_algorithm` | Identified by the general Pearl-Shpitser ID algorithm; finite-support plug-in estimation with EIF CIs is available for discrete variables |
| `:not_identified` | Not identified by the implemented criteria |

These strategies correspond to different identification arguments:

- `:a_fixable` covers backdoor-style effects. The treatment can be adjusted
  using its Markov pillow.
- `:p_fixable` covers front-door/primal-fixable/NPS effects. The treatment may
  be confounded with the outcome, but its children are not in the treatment's
  hidden-confounding district.
- `:nested_fixable` covers ADMGs where neither simple adjustment nor
  p-fixability is enough, but a One-Line ID style fixing sequence still
  identifies the effect.
- `:id_algorithm` covers more general ADMGs where the recursive ID algorithm
  finds a symbolic functional. The package can estimate such functionals with a
  finite-support plug-in estimator with EIF confidence intervals when the
  variables are discrete.

The returned strategy is an estimation instruction for `estimate_causal()` when
it is `:a_fixable`, `:p_fixable`, `:nested_fixable`, or `:id_algorithm`.
For `:id_algorithm`, inspect `identify(...).id_expression`, call
`ID_algorithm()` directly, or use `estimate_id()` for the finite-support
plug-in estimator.

## Optional NPCausal.jl Estimation

`CausalGraphs.jl` can also act as the identification layer for
[`NPCausal.jl`](https://github.com/xiangao/NPCausal.jl). Install and load both
packages, then call `estimate_causal_npcausal(...)`:

```julia
using Pkg
Pkg.add(url="https://github.com/xiangao/NPCausal.jl")

using CausalGraphs, NPCausal

res = estimate_causal_npcausal(
    a = [1, 0],
    data = df,
    graph = g_bd,
    treatment = :A,
    outcome = :Y,
    nsplits = 5,
)

r = x -> round(x, sigdigits=4)
(ACE = r(res[:NPCausalOnestep].ACE),
 SE = r(res[:NPCausalOnestep].standard_error))
```

The bridge currently supports `:a_fixable`/backdoor effects by passing the
treatment Markov pillow to `NPCausal.ate`. Front-door, nested-fixable, and
general ID-algorithm targets still use the estimators inside `CausalGraphs.jl`.

For the backdoor graph, the a-fixability condition holds because `A` is not in
the same hidden-confounding district as any descendant other than itself. The
Markov pillow is `X`, which becomes the adjustment set.

```@example graphs_identification
(descendants_A = descendants(g_bd, [:A]),
 district_A = district(g_bd, :A),
 markov_pillow_A = markov_pillow(g_bd, :A; treatment=:A))
```

If a graph has both `A -> Y` and `A <-> Y`, the directed causal effect and the
hidden common cause are entangled. That bow graph is the simplest
non-identified case.

## General ID Algorithm

`ID_algorithm()` runs the Pearl-Shpitser recursive ID algorithm for ADMG
queries. It is useful for separating mathematical identification from estimator
routing. On the front-door graph, it recovers the familiar front-door
functional:

```@example graphs_identification
id_general = ID_algorithm(g_fd, :A, :Y)
(identified = id_general.identified,
 expression = string(id_general.expression))
```

On the bow graph, the general ID algorithm fails and returns a hedge witness.

```@example graphs_identification
g_bow = make_graph(
    vertices = [:A, :Y],
    di_edges = [(:A, :Y)],
    bi_edges = [(:A, :Y)],
)

id_bow = ID_algorithm(g_bow, :A, :Y)
(identified = id_bow.identified, hedge = id_bow.hedge)
```
