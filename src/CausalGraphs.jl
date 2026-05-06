module CausalGraphs

using DataFrames
using GLM
using LinearAlgebra
using Random
using Statistics
using StatsModels
import GLM: predict

export ADMG, make_graph
export top_order, parents, children, descendants, ancestors
export district, all_districts, markov_blanket, markov_pillow
export subgraph, fixed_graph, reachable_closure
export is_fix, is_p_fix, is_np_saturated, mb_shielded
export to_dot, to_mermaid, to_svg, draw_graph, MermaidGraph, SVGGraph
export identify
export IDExpression, IDJoint, IDKernel, IDProduct, IDSum, ADMGIDResult
export id_algorithm, is_id_identified
export fixing_sequence, nested_fixability, is_nested_fixable
export superlearner
export estimate_causal, backdoor_tmle_a, nps_tmle_a, nested_anipw_a
export calculate_density_ratio_dnorm
export MDAG, MissingTree, IDLaw, IDResult, PropensityEstimate, CliqueResult, MEstimationResult
export make_mdag, make_mtree, append_leaf, prune_tree, are_trees_identical
export get_all_nodes, get_all_edges
export adj_matrix, top_order_R, nondescendants
export gfix, selectionx_set, problematic_set, colluder, colluder_chain, is_d_separated
export ID_algorithm, ID_flex_eval, summarize_ID
export propensity, find_clique, M_estimation, Mestimation
export f_ID_algorithm, f_ID_flex_eval, f_propensity, f_Mestimation
export r_indicator_for_x, x_variable_for_r
export compute_missing_weights

include("graph.jl")
include("missing_core.jl")
include("visualization.jl")
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
- general ID-algorithm identifiable: `identify()` reports `:id_algorithm`,
  but arbitrary ID functionals do not yet have an automatic estimator here

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
                            treatment=treatment, outcome=outcome,
                            sample_weights=sample_weights; kwargs...))
    elseif id.strategy == :p_fixable
        is_np_saturated(g) ||
            @info "Graph is not NP-saturated. More efficient estimators may exist."
        @info "Treatment is p-fixable (front-door/NPS). Using NPS-TMLE."
        return combine_levels(a, aval ->
            nps_tmle_a(a=aval, data=data, graph=g,
                       treatment=treatment, outcome=outcome,
                       sample_weights=sample_weights; kwargs...))
    elseif id.strategy == :nested_fixable
        @info "Treatment is nested-fixable. Using ANIPW."
        avals = collect(a isa AbstractVector ? a : [a])
        isempty(avals) && error("`a` must be a scalar or length-two vector.")
        length(avals) > 2 && error("Use scalar for E[Y(a)] or length-two for ACE.")
        out1 = nested_anipw_a(a=avals[1], data=data, graph=g,
                               treatment=treatment, outcome=outcome,
                               id_result=id, sample_weights=sample_weights; kwargs...)
        if length(avals) == 1
            return Dict{Symbol,Any}(
                :ANIPW    => (EYa=out1.ANIPW.estimated_psi,
                              lower_ci=out1.ANIPW.lower_ci,
                              upper_ci=out1.ANIPW.upper_ci,
                              EIF=out1.ANIPW.EIF),
                :NIPW     => (EYa=out1.NIPW.estimated_psi,),
                :ANIPW_Ya => out1.ANIPW,
                :rw       => out1.rw,
            )
        end
        out0 = nested_anipw_a(a=avals[2], data=data, graph=g,
                               treatment=treatment, outcome=outcome,
                               id_result=id, sample_weights=sample_weights; kwargs...)
        ace     = out1.ANIPW.estimated_psi - out0.ANIPW.estimated_psi
        eif_ace = out1.ANIPW.EIF - out0.ANIPW.EIF
        lo, hi  = ci_from_eif(ace, eif_ace, length(eif_ace))
        return Dict{Symbol,Any}(
            :ANIPW    => (ACE=ace, lower_ci=lo, upper_ci=hi, EIF=eif_ace),
            :NIPW     => (ACE=out1.NIPW.estimated_psi - out0.NIPW.estimated_psi,),
            :ANIPW_Y1 => out1.ANIPW,
            :ANIPW_Y0 => out0.ANIPW,
            :rw       => out1.rw,
        )
    else  # :id_algorithm
        error("The effect of $treatment on $outcome is identified by the general ID algorithm, " *
              "but `estimate_causal()` does not yet implement an estimator for arbitrary ID functionals. " *
              "Inspect `identify(...).id_expression` or `ID_algorithm(graph, treatment, outcome)`.")
    end
end

end # module
