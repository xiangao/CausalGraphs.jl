module CausalGraphsNPCausalExt

using CausalGraphs
using DataFrames
using NPCausal
using Statistics

function _npc_levels(a)
    avals = collect(a isa AbstractVector ? a : [a])
    isempty(avals) && error("`a` must be a scalar or length-two vector.")
    length(avals) > 2 && error("Use scalar for E[Y(a)] or length-two for ACE.")
    avals
end

function _npc_index(avals, a)
    idx = findfirst(x -> x == a, avals)
    idx === nothing &&
        error("Treatment level `$a` was not observed; NPCausal cannot estimate E[Y($a)].")
    idx
end

function _npc_mean_result(raw, idx::Int)
    est = raw.means.Estimate[idx]
    contrib = raw.ifvals[:, idx]
    eif = contrib .- est
    (EYa = est,
     lower_ci = raw.means.CI_Lower[idx],
     upper_ci = raw.means.CI_Upper[idx],
     standard_error = raw.means.StdError[idx],
     EIF = eif)
end

function _combine_npc_ate(a, raw, observed_avals)
    avals = _npc_levels(a)
    idx1 = _npc_index(observed_avals, avals[1])
    out1 = _npc_mean_result(raw, idx1)
    if length(avals) == 1
        return Dict{Symbol,Any}(
            :NPCausalOnestep => out1,
            :NPCausalRaw => raw,
        )
    end

    idx0 = _npc_index(observed_avals, avals[2])
    out0 = _npc_mean_result(raw, idx0)
    contrib = raw.ifvals[:, idx1] .- raw.ifvals[:, idx0]
    ace = mean(contrib)
    se = std(contrib) / sqrt(length(contrib))
    eif = contrib .- ace
    Dict{Symbol,Any}(
        :NPCausalOnestep => (ACE = ace,
                             lower_ci = ace - 1.96 * se,
                             upper_ci = ace + 1.96 * se,
                             standard_error = se,
                             EIF = eif),
        :NPCausalOnestep_Y1 => out1,
        :NPCausalOnestep_Y0 => out0,
        :NPCausalRaw => raw,
    )
end

function estimate_causal_npcausal_impl(; a, data::DataFrame,
                                         graph::Union{CausalGraphs.ADMG,Nothing}=nothing,
                                         vertices=nothing,
                                         di_edges=CausalGraphs.Edge[],
                                         bi_edges=CausalGraphs.Edge[],
                                         multivariate_variables=Dict{Symbol,Vector{Symbol}}(),
                                         treatment, outcome,
                                         sample_weights=nothing, kwargs...)
    sample_weights === nothing ||
        error("The NPCausal.jl bridge currently does not support `sample_weights`.")

    g = CausalGraphs.with_graph(graph, vertices, di_edges, bi_edges, multivariate_variables)
    id = CausalGraphs.identify(g, treatment, outcome)
    id.strategy == :a_fixable ||
        error("The NPCausal.jl bridge currently supports only `:a_fixable`/backdoor effects; got `$(id.strategy)`.")

    A = Symbol(treatment)
    Y = Symbol(outcome)
    predictors = CausalGraphs.replace_vector(
        CausalGraphs.markov_pillow(g, A; treatment=A),
        g.multivariate_variables,
    )
    all(p -> hasproperty(data, p), predictors) ||
        error("Data are missing at least one NPCausal adjustment variable.")

    X = isempty(predictors) ? DataFrame(_cg_intercept = ones(nrow(data))) : data[:, predictors]
    raw = NPCausal.ate(data[!, Y], data[!, A], X; kwargs...)
    observed_avals = sort(collect(unique(data[!, A])))
    _combine_npc_ate(a, raw, observed_avals)
end

end
