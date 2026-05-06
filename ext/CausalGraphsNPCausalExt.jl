module CausalGraphsNPCausalExt

using CausalGraphs
using DataFrames
using NPCausal

function estimate_causal_npcausal_impl(; a, data::DataFrame,
                                         graph::Union{CausalGraphs.ADMG,Nothing}=nothing,
                                         vertices=nothing,
                                         di_edges=CausalGraphs.Edge[],
                                         bi_edges=CausalGraphs.Edge[],
                                         multivariate_variables=Dict{Symbol,Vector{Symbol}}(),
                                         treatment, outcome,
                                         sample_weights=nothing, kwargs...)
    g = CausalGraphs.with_graph(graph, vertices, di_edges, bi_edges, multivariate_variables)
    NPCausal.admg_estimate_causal(a=a, data=data, graph=g,
                                  treatment=treatment, outcome=outcome,
                                  sample_weights=sample_weights; kwargs...)
end

end
