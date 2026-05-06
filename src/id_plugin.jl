# ── Generic finite-support plug-in estimator for symbolic ADMG ID ─────────────

function _id_supports(data::DataFrame, vars::Vector{Symbol}; max_levels::Int=25)
    supports = Dict{Symbol,Vector{Any}}()
    for v in vars
        hasproperty(data, v) || error("Data are missing graph variable `$v`.")
        any(ismissing, data[!, v]) &&
            error("`estimate_id` does not currently support missing values in graph variables.")
        vals = collect(unique(data[!, v]))
        length(vals) <= max_levels ||
            error("Variable `$v` has $(length(vals)) observed values. " *
                  "`estimate_id` is a finite-support plug-in estimator; " *
                  "increase `max_levels` only for genuinely discrete variables.")
        supports[v] = sort(vals)
    end
    supports
end

function _with_assignments(f, vars::Vector{Symbol}, supports, base::Dict{Symbol,Any})
    if isempty(vars)
        return f(base)
    end
    total = 0.0
    grids = (supports[v] for v in vars)
    for values in Iterators.product(grids...)
        a = copy(base)
        for (v, value) in zip(vars, values)
            a[v] = value
        end
        total += f(a)
    end
    total
end

function _empirical_joint_prob(data::DataFrame, vars::Vector{Symbol},
                               assignment::Dict{Symbol,Any}, weights::Vector{Float64})
    isempty(vars) && return 1.0
    for v in vars
        haskey(assignment, v) || error("No assigned value for `$v` while evaluating ID expression.")
    end
    num = 0.0
    for i in 1:nrow(data)
        ok = true
        for v in vars
            if data[!, v][i] != assignment[v]
                ok = false
                break
            end
        end
        ok && (num += weights[i])
    end
    num / sum(weights)
end

_expr_scope(expr::IDJoint) = copy(expr.vars)
_expr_scope(expr::IDKernel) = unique(vcat(expr.vars, expr.given))
_expr_scope(expr::IDProduct) = unique(vcat([_expr_scope(f) for f in expr.factors]...))
_expr_scope(expr::IDSum) = setdiff(_expr_scope(expr.expr), expr.vars)

function _restrict_assignment(assignment::Dict{Symbol,Any}, vars::Vector{Symbol})
    Dict{Symbol,Any}(v => assignment[v] for v in vars if haskey(assignment, v))
end

function _id_eval_source(expr::IDExpression, assignment, data, supports, weights)
    missing_vars = setdiff(_expr_scope(expr), collect(keys(assignment)))
    isempty(missing_vars) && return _id_eval(expr, assignment, data, supports, weights)
    _with_assignments(missing_vars, supports, copy(assignment)) do a
        _id_eval(expr, a, data, supports, weights)
    end
end

function _id_eval(expr::IDJoint, assignment, data, supports, weights)
    _empirical_joint_prob(data, expr.vars, assignment, weights)
end

function _id_eval(expr::IDKernel, assignment, data, supports, weights)
    numerator_assignment = _restrict_assignment(assignment, unique(vcat(expr.vars, expr.given)))
    numerator = _id_eval_source(expr.source, numerator_assignment, data, supports, weights)

    given_assignment = _restrict_assignment(assignment, expr.given)
    denominator = _with_assignments(expr.vars, supports, given_assignment) do a
        _id_eval_source(expr.source, a, data, supports, weights)
    end
    denominator <= eps() && return 0.0
    numerator / denominator
end

function _id_eval(expr::IDProduct, assignment, data, supports, weights)
    prod(_id_eval(f, assignment, data, supports, weights) for f in expr.factors)
end

function _id_eval(expr::IDSum, assignment, data, supports, weights)
    _with_assignments(expr.vars, supports, copy(assignment)) do a
        _id_eval(expr.expr, a, data, supports, weights)
    end
end

function _id_joint_cells(vars::Vector{Symbol}, supports)
    cells = Vector{Dict{Symbol,Any}}()
    for values in Iterators.product((supports[v] for v in vars)...)
        push!(cells, Dict{Symbol,Any}(v => value for (v, value) in zip(vars, values)))
    end
    cells
end

function _cell_matches(cell::Dict{Symbol,Any}, vars::Vector{Symbol}, assignment::Dict{Symbol,Any})
    for v in vars
        if cell[v] != assignment[v]
            return false
        end
    end
    true
end

function _id_cell_index(row::DataFrameRow, vars::Vector{Symbol}, cells)
    for (j, cell) in pairs(cells)
        all(v -> row[v] == cell[v], vars) && return j
    end
    error("Observed row contains a graph-variable combination outside the support.")
end

function _id_cell_probabilities(data::DataFrame, vars::Vector{Symbol}, cells, weights::Vector{Float64})
    p = zeros(Float64, length(cells))
    for i in 1:nrow(data)
        j = _id_cell_index(data[i, :], vars, cells)
        p[j] += weights[i]
    end
    p ./ sum(p)
end

function _joint_prob_from_cells(vars::Vector{Symbol}, assignment::Dict{Symbol,Any}, cells, p)
    isempty(vars) && return 1.0
    for v in vars
        haskey(assignment, v) || error("No assigned value for `$v` while evaluating ID expression.")
    end
    total = 0.0
    for (cell, prob) in zip(cells, p)
        _cell_matches(cell, vars, assignment) && (total += prob)
    end
    total
end

function _id_eval_source_p(expr::IDExpression, assignment, supports, cells, p)
    missing_vars = setdiff(_expr_scope(expr), collect(keys(assignment)))
    isempty(missing_vars) && return _id_eval_p(expr, assignment, supports, cells, p)
    _with_assignments(missing_vars, supports, copy(assignment)) do a
        _id_eval_p(expr, a, supports, cells, p)
    end
end

function _id_eval_p(expr::IDJoint, assignment, supports, cells, p)
    _joint_prob_from_cells(expr.vars, assignment, cells, p)
end

function _id_eval_p(expr::IDKernel, assignment, supports, cells, p)
    numerator_assignment = _restrict_assignment(assignment, unique(vcat(expr.vars, expr.given)))
    numerator = _id_eval_source_p(expr.source, numerator_assignment, supports, cells, p)

    given_assignment = _restrict_assignment(assignment, expr.given)
    denominator = _with_assignments(expr.vars, supports, given_assignment) do a
        _id_eval_source_p(expr.source, a, supports, cells, p)
    end
    denominator <= eps() && return 0.0
    numerator / denominator
end

function _id_eval_p(expr::IDProduct, assignment, supports, cells, p)
    prod(_id_eval_p(f, assignment, supports, cells, p) for f in expr.factors)
end

function _id_eval_p(expr::IDSum, assignment, supports, cells, p)
    _with_assignments(expr.vars, supports, copy(assignment)) do a
        _id_eval_p(expr.expr, a, supports, cells, p)
    end
end

function _id_psi_from_p(a, A::Symbol, Y::Symbol, expr::IDExpression, supports, cells, p)
    psi = 0.0
    probabilities = Dict{Any,Float64}()
    for y in supports[Y]
        assignment = Dict{Symbol,Any}(A => a, Y => y)
        py = _id_eval_p(expr, assignment, supports, cells, p)
        probabilities[y] = py
        psi += Float64(y) * py
    end
    psi, probabilities
end

function _id_effective_n(weights::Vector{Float64}, frequency_weights::Union{Bool,Nothing})
    use_frequency_weights = frequency_weights === nothing ?
        all(w -> isapprox(w, round(w); atol=1e-8), weights) :
        frequency_weights
    use_frequency_weights && return sum(weights)
    sum(weights)^2 / sum(weights .^ 2)
end

function _id_delta_eif(a, A::Symbol, Y::Symbol, expr::IDExpression, supports, cells, p;
                       eps_step::Float64=1e-6)
    psi, probabilities = _id_psi_from_p(a, A, Y, expr, supports, cells, p)
    eif = similar(p)
    for j in eachindex(p)
        p_eps = (1 - eps_step) .* p
        p_eps[j] += eps_step
        psi_eps, _ = _id_psi_from_p(a, A, Y, expr, supports, cells, p_eps)
        eif[j] = (psi_eps - psi) / eps_step
    end
    eif .-= sum(p .* eif)
    psi, probabilities, eif
end

function _id_eif_ci(estimate, cell_eif, p, n_eff)
    se = sqrt(sum(p .* (cell_eif .^ 2)) / n_eff)
    estimate - 1.96 * se, estimate + 1.96 * se, se
end

"""
    id_plugin_a(; a, data, graph, treatment, outcome, [id_result, max_levels])

Estimate `E[outcome(a)]` from a symbolic ADMG ID expression using a
finite-support plug-in estimator.

This estimator enumerates the observed support of the graph variables and
evaluates the ID functional using empirical conditional probabilities. It is
intended for discrete variables. It also computes a finite-support delta-method
EIF and Wald confidence interval for the plug-in functional. For continuous
variables, use one of the specialized semiparametric estimators when the graph is
a-, p-, or nested-fixable.

When `sample_weights` are supplied, integer-like weights are treated as
frequency weights for the IDPlugin standard error by default. Non-integer weights
use Kish effective sample size. Set `frequency_weights=true` or `false` to
override this behavior.
"""
function id_plugin_a(; a, data::DataFrame, graph::ADMG, treatment, outcome,
                       id_result=nothing, max_levels::Int=25,
                       sample_weights=nothing,
                       frequency_weights::Union{Bool,Nothing}=nothing,
                       eif_step::Float64=1e-6, kwargs...)
    A = sym(treatment)
    Y = sym(outcome)
    id = id_result === nothing ? ID_algorithm(graph, A, Y) : id_result
    id.identified || error("The effect of $treatment on $outcome is not identified by ID.")
    length(id.treatment) == 1 || error("`id_plugin_a` currently supports one treatment variable.")
    length(id.outcome) == 1 || error("`id_plugin_a` currently supports one outcome variable.")

    supports = _id_supports(data, graph.vertices; max_levels=max_levels)
    a in supports[A] ||
        error("Treatment level `$a` was not observed in `$A`; the finite-support plug-in estimator cannot evaluate it.")
    all(y -> y isa Number, supports[Y]) ||
        error("Outcome `$Y` must be numeric to estimate an expected value.")

    weights = observation_weights(sample_weights, nrow(data))
    cells = _id_joint_cells(graph.vertices, supports)
    p = _id_cell_probabilities(data, graph.vertices, cells, weights)
    estimated, probabilities, cell_eif =
        _id_delta_eif(a, A, Y, id.expression, supports, cells, p; eps_step=eif_step)
    total_probability = sum(values(probabilities))
    n_eff = _id_effective_n(weights, frequency_weights)
    lo, hi, se = _id_eif_ci(estimated, cell_eif, p, n_eff)
    row_cell_index = [_id_cell_index(data[i, :], graph.vertices, cells) for i in 1:nrow(data)]
    row_eif = cell_eif[row_cell_index]

    (estimated_psi = estimated,
     lower_ci = lo,
     upper_ci = hi,
     standard_error = se,
     EIF = row_eif,
     cell_EIF = cell_eif,
     cell_probabilities = p,
     cells = cells,
     n_eff = n_eff,
     probabilities = probabilities,
     total_probability = total_probability,
     expression = id.expression,
     id_result = id)
end

function _combine_id_levels(a, call_one)
    avals = collect(a isa AbstractVector ? a : [a])
    isempty(avals) && error("`a` must be a scalar or length-two vector.")
    length(avals) > 2 && error("Use scalar for E[Y(a)] or length-two for ACE.")
    out1 = call_one(avals[1])
    if length(avals) == 1
        return Dict{Symbol,Any}(
            :IDPlugin => (EYa=out1.estimated_psi,
                          lower_ci=out1.lower_ci,
                          upper_ci=out1.upper_ci,
                          standard_error=out1.standard_error,
                          EIF=out1.EIF,
                          total_probability=out1.total_probability,
                          expression=out1.expression),
            :IDPlugin_Ya => out1,
        )
    end
    out0 = call_one(avals[2])
    cell_eif = out1.cell_EIF - out0.cell_EIF
    row_eif = out1.EIF - out0.EIF
    ace = out1.estimated_psi - out0.estimated_psi
    lo, hi, se = _id_eif_ci(ace, cell_eif, out1.cell_probabilities, out1.n_eff)
    Dict{Symbol,Any}(
        :IDPlugin => (ACE=ace,
                      lower_ci=lo,
                      upper_ci=hi,
                      standard_error=se,
                      EIF=row_eif,
                      total_probability_a1=out1.total_probability,
                      total_probability_a0=out0.total_probability,
                      expression=out1.expression),
        :IDPlugin_Y1 => out1,
        :IDPlugin_Y0 => out0,
    )
end

"""
    estimate_id(; a, data, graph, treatment, outcome, [max_levels])

Run the general ADMG `ID_algorithm` and estimate the resulting symbolic
functional with the finite-support plug-in estimator.

The return value uses the `:IDPlugin` key. This is an automatic estimator for
discrete ID functionals. It includes a finite-support delta-method EIF and Wald
confidence interval for the plug-in functional. This is not a general TMLE/EIF
estimator for arbitrary continuous ADMG functionals.

When `sample_weights` are supplied, integer-like weights are treated as
frequency weights for the IDPlugin standard error by default. Non-integer weights
use Kish effective sample size. Set `frequency_weights=true` or `false` to
override this behavior.
"""
function estimate_id(; a, data::DataFrame, graph::ADMG, treatment, outcome,
                       sample_weights=nothing, kwargs...)
    id = ID_algorithm(graph, treatment, outcome)
    id.identified || error("The effect of $treatment on $outcome is not identified by ID.")
    _combine_id_levels(a, aval ->
        id_plugin_a(a=aval, data=data, graph=graph, treatment=treatment,
                    outcome=outcome, id_result=id,
                    sample_weights=sample_weights; kwargs...))
end
