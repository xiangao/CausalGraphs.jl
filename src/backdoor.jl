# ── Backdoor / a-fixable estimators ──────────────────────────────────────────

function backdoor_tmle_a(; a, data::DataFrame, graph::ADMG, treatment, outcome,
                           truncate_lower=0.0, truncate_upper=1.0,
                           formula_Y=nothing, formula_A=nothing,
                           superlearner_Y=false, superlearner_A=false,
                           crossfit=false, K=5, sample_weights=nothing, kwargs...)
    Aname, Yname = sym(treatment), sym(outcome)
    A, Y = col(data, Aname), col(data, Yname)
    n  = nrow(data)
    a0 = Float64(a); a1 = 1.0 - a0
    check_binary_treatment!(A, a0)
    sw = observation_weights(sample_weights, n)
    fit_sw = sample_weights === nothing ? nothing : sw

    mpA        = replace_vector(markov_pillow(graph, Aname; treatment=Aname),
                                graph.multivariate_variables)
    predictors = unique(vcat(mpA, [Aname]))
    mu_model   = superlearner_Y ? superlearner(binary=is_binary(Y)) : nothing
    pi_model   = superlearner_A ? superlearner(binary=true)         : nothing

    mu0, mu1 = if crossfit && mu_model !== nothing
        sample_weights === nothing ||
            error("`sample_weights` with cross-fitted SuperLearner outcome models is not currently supported.")
        crossfit_outcome(mu_model, data, Y, predictors, Aname, a0, a1, K)
    else
        fit_outcome_predictions(data, Y, predictors, Aname, a0, a1;
                                formula=formula_Y, response_name=Yname, ml_model=mu_model,
                                sample_weights=fit_sw)
    end
    mu_a = mu0

    pA1 = if crossfit && pi_model !== nothing
        sample_weights === nothing ||
            error("`sample_weights` with cross-fitted SuperLearner propensity models is not currently supported.")
        crossfit_propensity(pi_model, data, A, mpA, K)
    else
        fit_propensity(data, A, mpA; formula=formula_A, treatment_name=Aname,
                       ml_model=pi_model, sample_weights=fit_sw)
    end
    p_a = clip(probability_a(pA1, a0), truncate_lower, truncate_upper)

    # One-step / AIPW
    onestep_contrib = (A .== a0) .* (Y .- mu_a) ./ p_a .+ mu_a
    onestep_est = weighted_mean(onestep_contrib, sw)
    eif         = (sw ./ mean(sw)) .* (onestep_contrib .- onestep_est)
    lo, hi      = ci_from_eif(onestep_est, eif, n)
    onestep     = (estimated_psi=onestep_est, lower_ci=lo, upper_ci=hi,
                   EIF=eif, p_a_mpA=p_a, mu_Y_a=mu_a)

    gcomp = (estimated_psi=weighted_mean(mu_a, sw), mu_Y_a=mu_a)
    ipw_contrib = (A .== a0) .* Y ./ p_a
    ipw   = (estimated_psi=weighted_mean(ipw_contrib, sw), p_a_mpA=p_a)

    # TMLE targeting
    H = (A .== a0) ./ p_a
    if is_binary(Y)
        ε    = scalar_logistic_fluctuation(Y, safe_logit(mu_a), H; weights=sw)
        mu_t = expit(safe_logit(mu_a) .+ ε .* H)
    else
        ε    = weighted_mean_residual(Y, mu_a, sw .* H)
        mu_t = mu_a .+ ε
    end
    tmle_contrib = (A .== a0) .* (Y .- mu_t) ./ p_a .+ mu_t
    tmle_est = weighted_mean(tmle_contrib, sw)
    tmle_eif = (sw ./ mean(sw)) .* (tmle_contrib .- tmle_est)
    lo, hi   = ci_from_eif(tmle_est, tmle_eif, n)
    tmle     = (estimated_psi=tmle_est, lower_ci=lo, upper_ci=hi, EIF=tmle_eif,
                p_a_mpA=p_a, mu_Y_a=mu_t)

    return (TMLE=tmle, Onestep=onestep, Gcomp=gcomp, IPW=ipw)
end

# ── Shared: combine scalar/ACE results ───────────────────────────────────────
function combine_levels(a, call_one)
    avals = collect(a isa AbstractVector ? a : [a])
    isempty(avals) && error("`a` must be a scalar or length-two vector.")
    length(avals) > 2 && error("Use scalar for E[Y(a)] or length-two for ACE.")
    out1 = call_one(avals[1])
    if length(avals) == 1
        mk_eya(x) = (EYa=x.estimated_psi, lower_ci=x.lower_ci, upper_ci=x.upper_ci, EIF=x.EIF)
        out = Dict{Symbol,Any}(:TMLE    => mk_eya(out1.TMLE),
                                :Onestep => mk_eya(out1.Onestep),
                                :TMLE_Ya => out1.TMLE, :Onestep_Ya => out1.Onestep)
        hasproperty(out1, :Gcomp) && (out[:Gcomp] = (EYa=out1.Gcomp.estimated_psi,))
        hasproperty(out1, :IPW)   && (out[:IPW]   = (EYa=out1.IPW.estimated_psi,))
        return out
    end
    out0 = call_one(avals[2])
    function mk_ace(x1, x0)
        ate = x1.estimated_psi - x0.estimated_psi
        eif = x1.EIF - x0.EIF
        lo, hi = ci_from_eif(ate, eif, length(eif))
        (ACE=ate, lower_ci=lo, upper_ci=hi, EIF=eif)
    end
    out = Dict{Symbol,Any}(:TMLE    => mk_ace(out1.TMLE,    out0.TMLE),
                            :Onestep => mk_ace(out1.Onestep, out0.Onestep),
                            :TMLE_Y1 => out1.TMLE,  :TMLE_Y0 => out0.TMLE,
                            :Onestep_Y1 => out1.Onestep, :Onestep_Y0 => out0.Onestep)
    hasproperty(out1, :Gcomp) &&
        (out[:Gcomp] = (ACE=out1.Gcomp.estimated_psi - out0.Gcomp.estimated_psi,))
    hasproperty(out1, :IPW) &&
        (out[:IPW]   = (ACE=out1.IPW.estimated_psi   - out0.IPW.estimated_psi,))
    out
end
