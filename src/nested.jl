# ── Nested-fixable: ANIPW estimator ──────────────────────────────────────────
#
# Translated from anankeR's CausalEffect private methods:
#   .get_nested_rebalanced_weights, .fit_intrinsic_kernel, .augmented_nested_ipw
#
# Reference: Bhattacharya, Nabi & Shpitser (2022) "Semiparametric inference for
# causal effects in graphical models with hidden variables."

function fit_intrinsic_kernel(graph::ADMG, data::DataFrame, D::Vector{Symbol},
                               n_order::Vector{Symbol}, rebalance_prob::Vector{Float64})
    pa_D     = setdiff(parents(graph, D), D)
    anc_D    = ancestors(graph, vcat(D, pa_D))
    G        = subgraph(graph, anc_D)
    remaining = setdiff(G.vertices, vcat(D, pa_D))
    fixing_prob = copy(rebalance_prob)

    # Phase 1: iteratively fix vertices outside D ∪ pa_D until pa_D is a-fixable in G
    while !_set_fixable(G, pa_D).fixable
        fixed_this_pass = false
        fixable_V = nothing
        for V in remaining
            chs = children(G, V)
            if isempty(chs)
                # leaf node: fix without model (probability contribution = 1)
                G = fixed_graph(G, [V])
                filter!(!=(V), remaining)
                fixed_this_pass = true
                break
            elseif _is_fixable(G, V)
                fixable_V = V
            end
        end
        if !fixed_this_pass && fixable_V !== nothing
            V    = fixable_V
            mb_V = markov_blanket(G, V)
            w    = 1.0 ./ max.(fixing_prob, 1e-10)
            pV   = fit_density_prob(data, V, mb_V, w)
            fixing_prob = fixing_prob .* pV
            G = fixed_graph(G, [V])
            filter!(!=(V), remaining)
        elseif !fixed_this_pass
            # nothing fixable — proceed anyway to avoid infinite loop
            break
        end
    end

    # Phase 2: fit pa_D using Markov pillow under n_order
    pa_D_ordered = [v for v in n_order if v in pa_D]
    for V in pa_D_ordered
        mp_V = markov_pillow_with_order(graph, V, n_order)
        w    = 1.0 ./ max.(fixing_prob, 1e-10)
        pV   = fit_density_prob(data, V, mp_V, w)
        fixing_prob = fixing_prob .* pV
    end

    # Phase 3: kernel probability for variables in D
    kernel_prob = ones(nrow(data))
    D_ordered = [v for v in n_order if v in D]
    for V in D_ordered
        earlier_D = [d for d in D if findfirst(==(d), n_order) < findfirst(==(V), n_order)]
        ctx = vcat(pa_D, earlier_D)
        ctx = filter(!=(V), ctx)
        w   = 1.0 ./ max.(fixing_prob, 1e-10)
        pV  = fit_density_prob(data, V, ctx, w)
        kernel_prob = kernel_prob .* pV
    end

    kernel_prob
end

function get_nested_weights(graph::ADMG, data::DataFrame, Aname::Symbol,
                             modified_districts::Vector{Vector{Symbol}},
                             n_order::Vector{Symbol})
    rebalance_prob = ones(nrow(data))

    for D in modified_districts
        # Accumulate 1/p(V | mb(V, n_order)) for each V in D
        for V in D
            mp_V = markov_pillow_with_order(graph, V, n_order)
            pV   = fit_density_prob(data, V, mp_V, ones(nrow(data)))
            rebalance_prob = rebalance_prob ./ max.(pV, 1e-10)
        end
        # Multiply by intrinsic kernel φ_D
        kernel = fit_intrinsic_kernel(graph, data, D, n_order, 1.0 ./ rebalance_prob)
        rebalance_prob = rebalance_prob .* kernel
    end

    # rw = 1 / rebalance_prob (normalised so mean ≈ 1)
    rw = 1.0 ./ max.(rebalance_prob, 1e-10)
    rw ./ mean(rw)
end

function nested_anipw_a(; a, data::DataFrame, graph::ADMG, treatment, outcome, id_result,
                           truncate_lower=1e-3, truncate_upper=1-1e-3,
                           sample_weights=nothing, kwargs...)
    Aname = sym(treatment); Yname = sym(outcome)
    A = col(data, Aname); Y = col(data, Yname)
    n = nrow(data); a0 = Float64(a)
    check_binary_treatment!(A, a0)
    sw = observation_weights(sample_weights, n)

    ystar   = id_result.ystar
    Gystar  = id_result.Gystar
    n_order = id_result.n_order

    # Districts in G_Y* that intersect with district(T)
    dist_T = district(graph, Aname)
    modified_districts = [D for D in all_districts(Gystar)
                          if !isempty(intersect(D, dist_T))]

    # Compute nested rebalancing weights
    rw = get_nested_weights(graph, data, Aname, modified_districts, n_order)
    fit_weights = rw .* sw

    # Markov pillow of T under n_order
    mp_T = markov_pillow_with_order(graph, Aname, n_order)

    # Fit p(T | mp_T) with rw weights
    pA1 = if isempty(mp_T)
        w_sum = sum(fit_weights); fill(sum(fit_weights .* (A .== 1)) / w_sum, n)
    else
        fit_propensity(data, A, mp_T; sample_weights=fit_weights)
    end
    p_a = clip(probability_a(pA1, a0), truncate_lower, truncate_upper)

    # Fit E[Y | T, mp_T] with rw weights
    predictors = unique(vcat(mp_T, [Aname]))
    Yhat0, _   = fit_outcome_predictions(data, Y, predictors, Aname, a0, 1.0-a0;
                                          sample_weights=fit_weights)
    Yhat_a = Yhat0  # counterfactual prediction at T=a0

    # ANIPW doubly-robust formula (anankeR eq.)
    indicator  = Float64.(A .== a0)
    anipw_contrib = indicator ./ p_a .* (Y .- Yhat_a) .+ Yhat_a
    anipw_est  = weighted_mean(anipw_contrib, sw)
    eif        = (sw ./ mean(sw)) .* (anipw_contrib .- anipw_est)
    lo, hi     = ci_from_eif(anipw_est, eif, n)

    # IPW only (nested)
    nipw_est = weighted_mean(indicator ./ p_a .* Y, sw)

    return (ANIPW    = (estimated_psi=anipw_est, lower_ci=lo, upper_ci=hi, EIF=eif,
                        p_a_mpT=p_a, mu_Y_a=Yhat_a),
            NIPW     = (estimated_psi=nipw_est, p_a_mpT=p_a),
            rw       = rw)
end
