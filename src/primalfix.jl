# ── p-fixable / NPS-TMLE estimators ──────────────────────────────────────────

function calculate_density_ratio_dnorm(a0, M, graph::ADMG, treatment, data::DataFrame)
    a0 = Float64(a0); a1 = 1.0 - a0
    m = sym(M); Aname = sym(treatment)
    mp = replace_vector(markov_pillow(graph, m; treatment=Aname),
                        graph.multivariate_variables)
    Aname in mp || return ones(nrow(data))
    dat0 = set_column(data[:, mp], Aname, a0)
    dat1 = set_column(data[:, mp], Aname, a1)

    if haskey(graph.multivariate_variables, m)
        vars = graph.multivariate_variables[m]
        preds0 = Vector{Vector{Float64}}(); preds1 = Vector{Vector{Float64}}()
        residuals = zeros(nrow(data), length(vars))
        for (j,v) in enumerate(vars)
            y = col(data, v)
            is_binary(y) && error("dnorm ratio not supported for multivariate binary components.")
            β = fit_lm(data[:, mp], y, mp)
            push!(preds0, predict_lm(β, dat0, mp))
            push!(preds1, predict_lm(β, dat1, mp))
            residuals[:,j] = y - predict_lm(β, data[:, mp], mp)
        end
        Σ = cov(residuals)
        return [_mvnormal_pdf([data[i,v] for v in vars], [p[i] for p in preds0], Σ) /
                _mvnormal_pdf([data[i,v] for v in vars], [p[i] for p in preds1], Σ)
                for i in 1:nrow(data)]
    end

    y = col(data, m)
    if is_binary(y)
        β = fit_logistic(data[:, mp], y, mp)
        p0 = predict_logistic(β, dat0, mp)
        p1 = predict_logistic(β, dat1, mp)
        return [(y[i]==1 ? p0[i] : 1-p0[i]) / (y[i]==1 ? p1[i] : 1-p1[i]) for i in eachindex(y)]
    else
        β  = fit_lm(data[:, mp], y, mp)
        μ0 = predict_lm(β, dat0, mp)
        μ1 = predict_lm(β, dat1, mp)
        σ  = std(y - predict_lm(β, data[:, mp], mp))
        return normal_pdf.(y, μ0, σ) ./ normal_pdf.(y, μ1, σ)
    end
end

function _mvnormal_pdf(x::Vector{Float64}, μ::Vector{Float64}, Σ::Matrix{Float64})
    k = length(x); S = Symmetric(Σ + 1e-8I); d = x - μ
    exp(-0.5*dot(d, S\d)) / sqrt((2π)^k * max(det(S), 1e-12))
end

function bayes_density_ratio(v::Symbol, a0, graph::ADMG, treatment::Symbol,
                              data::DataFrame, densratio_A;
                              zerodiv_avoid=0.0, ml_model=nothing)
    A = col(data, treatment)
    mp = markov_pillow(graph, v; treatment=treatment)
    cols_full = setdiff(replace_vector(unique(vcat([v], mp)),
                                       graph.multivariate_variables), [treatment])
    cols_mp   = setdiff(replace_vector(mp, graph.multivariate_variables), [treatment])
    p1_full   = fit_propensity(data, A, cols_full; ml_model)
    p0_full   = probability_a(p1_full, Float64(a0))
    p1a_full  = 1 .- p0_full
    p0_full   = max.(p0_full, zerodiv_avoid); p1a_full = max.(p1a_full, zerodiv_avoid)
    numerator = p0_full ./ p1a_full

    mpA = replace_vector(markov_pillow(graph, treatment; treatment=treatment),
                         graph.multivariate_variables)
    if same_set(cols_mp, mpA)
        return (ratio=numerator ./ densratio_A, tied_to_A=true, raw=numerator)
    end
    p1_mp  = fit_propensity(data, A, cols_mp; ml_model)
    p0_mp  = probability_a(p1_mp, Float64(a0))
    p1a_mp = 1 .- p0_mp
    p0_mp  = max.(p0_mp, zerodiv_avoid); p1a_mp = max.(p1a_mp, zerodiv_avoid)
    (ratio=numerator ./ (p0_mp ./ p1a_mp), tied_to_A=false, raw=numerator)
end

function fit_seq_reg!(mu::Dict{Tuple{Symbol,Int},Vector{Float64}}, v::Symbol,
                       tau::Vector{Symbol}, L::Vector{Symbol}, graph::ADMG,
                       treatment::Symbol, data::DataFrame, Y_is_binary::Bool,
                       a0, a1; ml_model=nothing)
    mpv  = replace_vector(markov_pillow(graph, v; treatment=treatment),
                          graph.multivariate_variables)
    nextv   = tau[findfirst(==(v), tau) + 1]
    nextmu  = nextv in L ? mu[(nextv,1)] : mu[(nextv,0)]
    response = Y_is_binary ? safe_logit(nextmu) : nextmu
    dat  = data[:, mpv]
    dat0 = treatment in mpv ? set_column(dat, treatment, a0) : dat
    dat1 = treatment in mpv ? set_column(dat, treatment, a1) : dat
    if ml_model !== nothing
        mach = _mlj_fit_regressor(ml_model, dat, response)
        reg0 = _mlj_predict_reg(mach, dat0); reg1 = _mlj_predict_reg(mach, dat1)
    else
        β    = fit_lm(dat, response, mpv)
        reg0 = predict_lm(β, dat0, mpv); reg1 = predict_lm(β, dat1, mpv)
    end
    mu[(v,0)] = Y_is_binary ? expit(reg0) : reg0
    mu[(v,1)] = Y_is_binary ? expit(reg1) : reg1
    nothing
end

function nps_tmle_a(; a, data::DataFrame, graph::ADMG, treatment, outcome,
                      ratio_method_L="bayes", ratio_method_M="bayes",
                      n_iter=500, cvg_criteria=0.01,
                      truncate_lower=0.0, truncate_upper=1.0,
                      zerodiv_avoid=0.0, formula_Y=nothing, formula_A=nothing,
                      superlearner_Y=false, superlearner_A=false, superlearner_seq=false,
                      superlearner_M=false, superlearner_L=false,
                      crossfit=false, K=5, kwargs...)
    (ratio_method_L == "densratio" || ratio_method_M == "densratio") &&
        error("Kernel densratio not implemented.")

    Aname, Yname = sym(treatment), sym(outcome)
    A, Y = col(data, Aname), col(data, Yname)
    n = nrow(data); a0 = Float64(a); a1 = 1.0 - a0
    tau   = top_order(graph; treatment=Aname)
    order = Dict(v => i for (i,v) in enumerate(tau))
    sets  = cml(graph, Aname)
    C     = rerank(sets.C, tau); L = rerank(sets.L, tau)
    L_noA = [v for v in L if v != Aname]; M = rerank(sets.M, tau)
    Y_is_binary = is_binary(Y)
    mu_model  = superlearner_Y   ? superlearner(binary=Y_is_binary) : nothing
    pi_model  = superlearner_A   ? superlearner(binary=true)        : nothing
    seq_model = superlearner_seq ? superlearner(binary=false)        : nothing
    M_model   = superlearner_M   ? superlearner(binary=true)         : nothing
    L_model   = superlearner_L   ? superlearner(binary=true)         : nothing

    mu  = Dict{Tuple{Symbol,Int},Vector{Float64}}()
    mpY = replace_vector(markov_pillow(graph, Yname; treatment=Aname),
                         graph.multivariate_variables)
    mu0, mu1 = if crossfit && mu_model !== nothing
        crossfit_outcome(mu_model, data, Y, mpY, Aname, a0, a1, K)
    else
        fit_outcome_predictions(data, Y, mpY, Aname, a0, a1;
                                formula=formula_Y, response_name=Yname, ml_model=mu_model)
    end
    mu[(Yname,0)] = mu0; mu[(Yname,1)] = mu1

    mpA = replace_vector(markov_pillow(graph, Aname; treatment=Aname),
                         graph.multivariate_variables)
    _raw_pA1 = if crossfit && pi_model !== nothing
        crossfit_propensity(pi_model, data, A, mpA, K)
    else
        fit_propensity(data, A, mpA; formula=formula_A, treatment_name=Aname, ml_model=pi_model)
    end
    pA1  = clip(_raw_pA1, truncate_lower, truncate_upper)
    p_a1 = probability_a(pA1, a1); p_a0 = 1 .- p_a1
    p_a0 = max.(p_a0, zerodiv_avoid); p_a1 = max.(p_a1, zerodiv_avoid)

    dens = Dict{Symbol,Vector{Float64}}(Aname => p_a0 ./ p_a1)
    raw_for_A_update = Dict{Symbol,Vector{Float64}}()
    tied_L = Symbol[]; tied_M = Symbol[]

    for v in L_noA
        if ratio_method_L == "dnorm"
            dens[v] = calculate_density_ratio_dnorm(a0, v, graph, Aname, data)
        elseif ratio_method_L == "bayes"
            r = bayes_density_ratio(v, a0, graph, Aname, data, dens[Aname];
                                    zerodiv_avoid, ml_model=L_model)
            dens[v] = r.ratio
            if r.tied_to_A; raw_for_A_update[v] = r.raw; push!(tied_L, v); end
        else; error("Invalid ratio_method_L."); end
    end

    M_includeA = [v for v in M if Aname in markov_pillow(graph, v; treatment=Aname)]
    for v in setdiff(M, M_includeA); dens[v] = ones(n); end
    for v in M_includeA
        if ratio_method_M == "dnorm"
            dens[v] = calculate_density_ratio_dnorm(a0, v, graph, Aname, data)
        elseif ratio_method_M == "bayes"
            r = bayes_density_ratio(v, a0, graph, Aname, data, dens[Aname];
                                    zerodiv_avoid, ml_model=M_model)
            dens[v] = r.ratio
            if r.tied_to_A; raw_for_A_update[v] = r.raw; push!(tied_M, v); end
        else; error("Invalid ratio_method_M."); end
    end

    between = order[Aname]+1 <= order[Yname]-1 ? tau[order[Aname]+1:order[Yname]-1] : Symbol[]
    for v in reverse(between)
        fit_seq_reg!(mu, v, tau, L, graph, Aname, data, Y_is_binary, a0, a1; ml_model=seq_model)
    end

    function compute_eifs!()
        sel_M = [m for m in M if order[m] < order[Yname]]
        sel_L = [l for l in L if order[l] < order[Yname]]
        fM = product_columns(dens, sel_M, n); fL = product_columns(dens, sel_L, n)
        EIFY = Yname in L ? (A .== a1) .* fM .* (Y .- mu[(Yname,1)]) :
                             (A .== a0) .* (1 ./ fL) .* (Y .- mu[(Yname,0)])
        EIFv = Dict{Symbol,Vector{Float64}}()
        for v in reverse(between)
            nextv  = tau[order[v]+1]
            nextmu = nextv in L ? mu[(nextv,1)] : mu[(nextv,0)]
            if v in L
                f = product_columns(dens, [m for m in M if order[m] < order[v]], n)
                EIFv[v] = (A .== a1) .* f .* (nextmu .- mu[(v,1)])
            else
                f = product_columns(dens, [l for l in L if order[l] < order[v]], n)
                EIFv[v] = (A .== a0) .* (1 ./ f) .* (nextmu .- mu[(v,0)])
            end
        end
        nextA = tau[order[Aname]+1]
        EIFA  = ((A .== a1) .- p_a1) .* mu[(nextA,0)]
        sumv  = isempty(between) ? zeros(n) : reduce(+, [EIFv[v] for v in between])
        EIFY, EIFA, sumv, nextA
    end

    EIFY, EIFA, sumv, nextA = compute_eifs!()
    estimate = mean(EIFY .+ sumv .+ EIFA .+ p_a1 .* mu[(nextA,0)] .+ (A .== a0) .* Y)
    eif = EIFY .+ sumv .+ EIFA .+ p_a1 .* mu[(nextA,0)] .+ (A .== a0) .* Y .- estimate
    lo, hi = ci_from_eif(estimate, eif, n)
    onestep = (estimated_psi=estimate, lower_ci=lo, upper_ci=hi, EIF=eif)

    EDstar = 10.0; EDrecord = Float64[]; iter = 0
    while abs(EDstar) > cvg_criteria && iter < n_iter
        cleverA = mu[(nextA,0)]
        εA  = scalar_logistic_fluctuation(A .== a1, safe_logit(p_a1), cleverA)
        p_a1 = clip(expit(safe_logit(p_a1) .+ εA .* cleverA), 1e-8, 1-1e-8)
        p_a0 = max.(1 .- p_a1, zerodiv_avoid); p_a1 = max.(p_a1, zerodiv_avoid)
        dens[Aname] = p_a0 ./ p_a1
        for v in tied_L; dens[v] = raw_for_A_update[v] ./ dens[Aname]; end
        for v in tied_M; dens[v] = raw_for_A_update[v] ./ dens[Aname]; end

        sel_M = [m for m in M if order[m] < order[Yname]]
        sel_L = [l for l in L if order[l] < order[Yname]]
        fM = product_columns(dens, sel_M, n); fL = product_columns(dens, sel_L, n)
        weightY = Yname in L ? (A .== a1) .* fM : (A .== a0) .* (1 ./ fL)
        offsetY = Yname in L ? mu[(Yname,1)] : mu[(Yname,0)]
        key_Y   = Yname in L ? (Yname,1) : (Yname,0)
        if Y_is_binary
            εY = scalar_logistic_fluctuation(Y, safe_logit(offsetY), weightY)
            mu[key_Y] = expit(safe_logit(offsetY) .+ εY .* weightY)
        else
            εY = weighted_mean_residual(Y, offsetY, weightY)
            mu[key_Y] = offsetY .+ εY
        end

        for v in reverse(between)
            fit_seq_reg!(mu, v, tau, L, graph, Aname, data, Y_is_binary, a0, a1; ml_model=seq_model)
            sel_Mv = [m for m in M if order[m] < order[v]]
            sel_Lv = [l for l in L if order[l] < order[v]]
            weightv = v in L ? (A .== a1) .* product_columns(dens, sel_Mv, n) :
                                (A .== a0) .* (1 ./ product_columns(dens, sel_Lv, n))
            nextv   = tau[order[v]+1]
            nextmu  = nextv in L ? mu[(nextv,1)] : mu[(nextv,0)]
            key     = v in L ? (v,1) : (v,0)
            offsetv = mu[key]
            if Y_is_binary
                εv = scalar_logistic_fluctuation(nextmu, safe_logit(offsetv), weightv)
                mu[key] = expit(safe_logit(offsetv) .+ εv .* weightv)
            else
                εv = weighted_mean_residual(nextmu, offsetv, weightv)
                mu[key] = offsetv .+ εv
            end
        end

        EIFY, EIFA, sumv, nextA = compute_eifs!()
        EDstar = mean(EIFA) + mean(EIFY) + mean(sumv)
        push!(EDrecord, EDstar); iter += 1
    end

    estimate = mean(EIFY .+ sumv .+ EIFA .+ p_a1 .* mu[(nextA,0)] .+ (A .== a0) .* Y)
    eif = EIFY .+ sumv .+ EIFA .+ p_a1 .* mu[(nextA,0)] .+ (A .== a0) .* Y .- estimate
    lo, hi = ci_from_eif(estimate, eif, n)
    tmle = (estimated_psi=estimate, lower_ci=lo, upper_ci=hi, EIF=eif,
            EDstar=EDstar, iter=iter, EDstar_record=EDrecord)
    (Onestep=onestep, TMLE=tmle)
end
