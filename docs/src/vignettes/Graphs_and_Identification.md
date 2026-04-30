# Graphs and Identification

```@meta
CurrentModule = CausalGraphs
```

This vignette introduces the graph layer and the `identify()` function.

## Constructing an ADMG

An acyclic directed mixed graph has directed edges for observed causal
relationships and bidirected edges for latent common causes.

```@example graphs_identification
using CausalGraphs

g_fd = make_graph(
    vertices = [:A, :M, :Y],
    di_edges = [(:A, :M), (:M, :Y)],
    bi_edges = [(:A, :Y)],
)

to_mermaid(g_fd; direction="LR")
```

```@example graphs_identification
g_bd = make_graph(
    vertices = [:X, :A, :Y],
    di_edges = [(:X, :A), (:X, :Y), (:A, :Y)],
)

to_mermaid(g_bd; direction="LR")
```

## Graph Properties

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

