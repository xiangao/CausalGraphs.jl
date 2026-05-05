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

The returned strategy is both an identification statement and an estimation
instruction for `estimate_causal()`.
