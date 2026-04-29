module CausalGraphs

using DataFrames
using GLM
using LinearAlgebra
using Random
using Statistics
using StatsModels

export ADMG, make_graph
export top_order, parents, children, descendants, ancestors
export district, all_districts, markov_blanket, markov_pillow
export subgraph, fixed_graph, reachable_closure
export is_fix, is_p_fix, is_np_saturated, mb_shielded
export identify
export superlearner
export estimate_causal, backdoor_tmle_a, nps_tmle_a, nested_anipw_a
export calculate_density_ratio_dnorm
export compute_missing_weights

include("graph.jl")
include("identification.jl")
include("nuisance.jl")
include("backdoor.jl")
include("primalfix.jl")
include("nested.jl")
include("missing.jl")
include("ensemble.jl")

# ── Unified entry point ───────────────────────────────────────────────────────

"""
    estimate_causal(; a, data, graph=nothing, vertices, di_edges, bi_edges,
                     treatment, outcome, [superlearner_Y, superlearner_A,
                     crossfit, K, sample_weights, kwargs...])

Auto-routing causal effect estimator for ADMGs.

Identification is determined automatically:
- a-fixable (backdoor): returns TMLE, Onestep, Gcomp, IPW
- p-fixable (front-door / NPS): returns TMLE, Onestep
- nested-fixable: returns ANIPW, NIPW

`a` may be a scalar for E[Y(a)] or a length-2 vector `[a1, a0]` for the ACE.

# Example — backdoor (a-fixable)
```julia
g = make_graph(vertices=[:A,:Y,:X],
               di_edges=[(:X,:A),(:X,:Y),(:A,:Y)])
res = estimate_causal(a=[1,0], data=df, graph=g, treatment=:A, outcome=:Y)
res[:TMLE].ACE
```

# Example — front-door (p-fixable)
```julia
g = make_graph(vertices=[:A,:M,:Y],
               bi_edges=[(:A,:Y)],
               di_edges=[(:A,:M),(:M,:Y)])
res = estimate_causal(a=[1,0], data=df, graph=g, treatment=:A, outcome=:Y)
```

# Example — nested-fixable (with SuperLearner)
```julia
res = estimate_causal(a=[1,0], data=df, graph=g, treatment=:A, outcome=:Y,
                      superlearner_Y=true, superlearner_A=true)
```
"""
function estimate_causal(; a, data::DataFrame,
                           graph::Union{ADMG,Nothing}=nothing,
                           vertices=nothing, di_edges=Edge[], bi_edges=Edge[],
                           multivariate_variables=Dict{Symbol,Vector{Symbol}}(),
                           treatment, outcome,
                           sample_weights=nothing, kwargs...)
    g  = with_graph(graph, vertices, di_edges, bi_edges, multivariate_variables)
    id = identify(g, treatment, outcome)

    if id.strategy == :not_identified
        error("The effect of $treatment on $outcome is not identified in this ADMG.")
    elseif id.strategy == :a_fixable
        @info "Treatment is a-fixable (backdoor). Using TMLE."
        return combine_levels(a, aval ->
            backdoor_tmle_a(a=aval, data=data, graph=g,
                            treatment=treatment, outcome=outcome; kwargs...))
    elseif id.strategy == :p_fixable
        is_np_saturated(g) ||
            @info "Graph is not NP-saturated. More efficient estimators may exist."
        @info "Treatment is p-fixable (front-door/NPS). Using NPS-TMLE."
        return combine_levels(a, aval ->
            nps_tmle_a(a=aval, data=data, graph=g,
                       treatment=treatment, outcome=outcome; kwargs...))
    else  # :nested_fixable
        @info "Treatment is nested-fixable. Using ANIPW."
        return combine_levels(a, aval ->
            nested_anipw_a(a=aval, data=data, graph=g,
                           treatment=treatment, outcome=outcome,
                           id_result=id; kwargs...))
    end
end

end # module
