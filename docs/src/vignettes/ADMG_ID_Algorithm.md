# General ADMG ID Algorithm

```@meta
CurrentModule = CausalGraphs
```

This vignette introduces the general ADMG identification tools in
`CausalGraphs.jl`. `identify()` is the estimator-routing function: it checks
the effect classes for which `estimate_causal()` has estimators. `ID_algorithm()`
is the more general symbolic Pearl-Shpitser ID check for ADMG total-effect
queries.

## Front-Door as ID

The front-door graph has hidden treatment-outcome confounding, represented by
`A <-> Y`, but the effect is still identified through the mediator `M`.

```@example admg_id
using CausalGraphs

g_fd = make_graph(
    vertices = [:A, :M, :Y],
    di_edges = [(:A, :M), (:M, :Y)],
    bi_edges = [(:A, :Y)],
)

draw_graph(g_fd; direction="LR")
```

For estimation, this graph is p-fixable, so `identify()` routes to the NPS-TMLE
estimator.

```@example admg_id
identify(g_fd, :A, :Y).strategy
```

The same graph is also identified by the general ID algorithm:

```@example admg_id
id_fd = ID_algorithm(g_fd, :A, :Y)
(identified = id_fd.identified,
 expression = string(id_fd.expression))
```

The printed expression is intentionally compact. In the front-door graph it is
the familiar formula: average the mediator distribution under the intervention,
then average the mediator-outcome regression over the observed treatment
distribution.

## Non-Identification

The bow graph has both `A -> Y` and hidden confounding `A <-> Y`. The observed
association cannot separate the directed effect from the latent common cause.

```@example admg_id
g_bow = make_graph(
    vertices = [:A, :Y],
    di_edges = [(:A, :Y)],
    bi_edges = [(:A, :Y)],
)

draw_graph(g_bow; direction="LR")
```

```@example admg_id
id_bow = ID_algorithm(g_bow, :A, :Y)
(identified = id_bow.identified,
 hedge = id_bow.hedge)
```

When ID fails, the returned hedge is a graph-theoretic witness of
non-identification. In practice, this is more useful than a bare `false`
because it points to the district structure responsible for the failure.

## Fixing Sequences

Fixing is the graph operation that turns random vertices into intervention-like
or conditioned-on vertices by removing incoming directed edges and bidirected
edges involving the fixed node. A vertex is fixable when it is not in the same
district as one of its proper descendants.

```@example admg_id
g_bd = make_graph(
    vertices = [:X, :A, :Y],
    di_edges = [(:X, :A), (:X, :Y), (:A, :Y)],
)

fixing_sequence(g_bd, [:A])
```

For the backdoor graph, fixing `A` succeeds immediately because `A` is not
bidirected-connected to its descendant `Y`.

## Nested-Fixability Diagnostics

Nested-fixability uses fixing operations inside induced subgraphs. It is more
general than ordinary backdoor adjustment and p-fixability, and it underlies the
ANIPW/NIPW estimator route in this package.

```@example admg_id
g_nested = make_graph(
    vertices = [:Z, :A, :M, :Y],
    di_edges = [(:Z, :A), (:Z, :Y), (:A, :M), (:M, :Y)],
    bi_edges = [(:A, :M)],
)

nested = nested_fixability(g_nested, :A, :Y)
(p_fixable = is_p_fix(g_nested, :A),
 nested_identified = nested.identified,
 ystar = nested.ystar,
 nested_order = nested.n_order)
```

The nested diagnostics show which outcome ancestors remain after fixing
treatment and which districts must be reachable by fixing vertices outside the
district.

## Practical Boundary

Use the functions at different levels of generality:

- `is_fix()` and `is_p_fix()` are simple graph checks.
- `nested_fixability()` explains the nested/fixing route used by the ANIPW
  estimator.
- `ID_algorithm()` answers the more general symbolic ADMG identification
  question and returns a hedge when the effect is not identified.
- `identify()` combines these checks for package routing.
- `estimate_causal()` estimates only the currently implemented estimator
  classes: a-fixable, p-fixable, and nested-fixable.

General ID can produce complicated functionals. Those functionals are useful
identification results, but they require separate estimation machinery before
they should be treated as automatic estimators.
