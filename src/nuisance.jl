# ── Nuisance fitting helpers ──────────────────────────────────────────────────
# Shared by backdoor, primalfix, and nested estimators.

function fit_outcome_predictions(data::DataFrame, y, predictors::Vector{Symbol},
                                  intervention_col::Symbol, a0, a1;
                                  formula=nothing, response_name::Union{Symbol,Nothing}=nothing,
                                  ml_model=nothing, sample_weights=nothing)
    dat  = data[:, predictors]
    dat0 = intervention_col in predictors ? set_column(dat, intervention_col, a0) : dat
    dat1 = intervention_col in predictors ? set_column(dat, intervention_col, a1) : dat
    if ml_model !== nothing
        if is_binary(y)
            mach = _mlj_fit_classifier(ml_model, dat, y)
            return (_mlj_predict_prob1(mach, dat0), _mlj_predict_prob1(mach, dat1))
        else
            mach = _mlj_fit_regressor(ml_model, dat, y)
            return (_mlj_predict_reg(mach, dat0), _mlj_predict_reg(mach, dat1))
        end
    end
    if formula !== nothing
        using_glm = true
        yname = response_name !== nothing ? response_name :
                (formula.lhs isa Term ? formula.lhs.sym :
                 error("Pass response_name explicitly."))
        dat_fit = copy(dat); dat_fit[!, yname] = Float64.(y)
        family = is_binary(y) ? Binomial() : Normal()
        link   = is_binary(y) ? LogitLink() : IdentityLink()
        fit_glm = glm(formula, dat_fit, family, link)
        return (Float64.(predict(fit_glm, dat0)), Float64.(predict(fit_glm, dat1)))
    end
    if is_binary(y)
        β = fit_logistic(dat, y, predictors; sample_weights)
        return (predict_logistic(β, dat0, predictors), predict_logistic(β, dat1, predictors))
    else
        β = fit_lm(dat, y, predictors; weights=sample_weights)
        return (predict_lm(β, dat0, predictors), predict_lm(β, dat1, predictors))
    end
end

function fit_propensity(data::DataFrame, A, predictors::Vector{Symbol};
                         formula=nothing, treatment_name::Union{Symbol,Nothing}=nothing,
                         ml_model=nothing, sample_weights=nothing)
    isempty(predictors) && return fill(mean(Float64.(A .== 1)), nrow(data))
    if ml_model !== nothing
        mach = _mlj_fit_classifier(ml_model, data[:, predictors], A)
        return _mlj_predict_prob1(mach, data[:, predictors])
    end
    if formula !== nothing
        aname = treatment_name !== nothing ? treatment_name :
                (formula.lhs isa Term ? formula.lhs.sym : error("Pass treatment_name."))
        dat = copy(data[:, predictors]); dat[!, aname] = Float64.(A)
        fit_glm = glm(formula, dat, Binomial(), LogitLink())
        return clip(Float64.(predict(fit_glm, data[:, predictors])), 1e-8, 1-1e-8)
    end
    β = fit_logistic(data[:, predictors], A, predictors; sample_weights)
    predict_logistic(β, data[:, predictors], predictors)
end

# ── Density estimation for nested case ───────────────────────────────────────
# Returns p(V = observed_value | predictors) per row, fitted with sample weights.

function fit_density_prob(data::DataFrame, V::Symbol, predictors::Vector{Symbol},
                           weights::Vector{Float64})
    y = col(data, V)
    if is_binary(y)
        if isempty(predictors)
            p_hat = fill(sum(y .* weights) / sum(weights), nrow(data))
        else
            β = fit_logistic(data[:, predictors], y, predictors; sample_weights=weights)
            p_hat = predict_logistic(β, data[:, predictors], predictors)
        end
        return y .* p_hat .+ (1 .- y) .* (1 .- p_hat)
    else
        if isempty(predictors)
            μ = fill(sum(y .* weights) / sum(weights), nrow(data))
        else
            β = fit_lm(data[:, predictors], y, predictors; weights)
            μ = predict_lm(β, data[:, predictors], predictors)
        end
        resid = y - μ
        w_sum = sum(weights)
        σ = sqrt(max(sum(weights .* resid .^ 2) / w_sum, 1e-8))
        return normal_pdf.(y, μ, σ)
    end
end

# ── Cross-fitting helpers (for a- and p-fixable with SuperLearner) ────────────
function _kfold_indices(n::Int, K::Int)
    idx = shuffle(1:n)
    [idx[round(Int,(k-1)*n/K)+1 : round(Int,k*n/K)] for k in 1:K]
end

function crossfit_propensity(model, data::DataFrame, A, predictors::Vector{Symbol}, K::Int)
    n = nrow(data); pA1 = zeros(n)
    folds = _kfold_indices(n, K)
    for vf in 1:K
        test  = folds[vf]
        train = vcat(folds[setdiff(1:K, vf)]...)
        mach  = _mlj_fit_classifier(model, data[train, predictors], A[train])
        pA1[test] = _mlj_predict_prob1(mach, data[test, predictors])
    end
    clip(pA1, 1e-8, 1-1e-8)
end

function crossfit_outcome(model, data::DataFrame, y, predictors::Vector{Symbol},
                           intervention_col::Symbol, a0, a1, K::Int)
    n = nrow(data); mu0 = zeros(n); mu1 = zeros(n)
    dat     = data[:, predictors]
    dat0all = intervention_col in predictors ? set_column(dat, intervention_col, a0) : dat
    dat1all = intervention_col in predictors ? set_column(dat, intervention_col, a1) : dat
    folds   = _kfold_indices(n, K)
    binary  = is_binary(y)
    for vf in 1:K
        test  = folds[vf]
        train = vcat(folds[setdiff(1:K, vf)]...)
        if binary
            mach = _mlj_fit_classifier(model, dat[train,:], y[train])
            mu0[test] = _mlj_predict_prob1(mach, dat0all[test,:])
            mu1[test] = _mlj_predict_prob1(mach, dat1all[test,:])
        else
            mach = _mlj_fit_regressor(model, dat[train,:], y[train])
            mu0[test] = _mlj_predict_reg(mach, dat0all[test,:])
            mu1[test] = _mlj_predict_reg(mach, dat1all[test,:])
        end
    end
    (mu0, mu1)
end
