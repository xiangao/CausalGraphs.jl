# ── Pearl-Shpitser ID algorithm for ADMGs ─────────────────────────────────────
#
# This is a symbolic implementation of the recursive ID algorithm for
# interventional distributions in semi-Markovian ADMGs. It reports whether
# P(outcome | do(treatment)) is identifiable and returns a readable symbolic
# functional. Estimation for arbitrary ID functionals is not routed through
# `estimate_causal()` yet.

abstract type IDExpression end

struct IDJoint <: IDExpression
    vars::Vector{Symbol}
end

struct IDKernel <: IDExpression
    vars::Vector{Symbol}
    given::Vector{Symbol}
    source::IDExpression
end

struct IDProduct <: IDExpression
    factors::Vector{IDExpression}
end

struct IDSum <: IDExpression
    vars::Vector{Symbol}
    expr::IDExpression
end

struct ADMGIDResult
    identified::Bool
    treatment::Vector{Symbol}
    outcome::Vector{Symbol}
    expression::Union{IDExpression,Nothing}
    hedge::Union{NamedTuple,Nothing}
end

struct _IDFailure
    hedge::NamedTuple
end

_ordered_intersect(xs, ys) = [x for x in symvec(xs) if x in Set(symvec(ys))]
_ordered_union(xs, ys) = unique(vcat(symvec(xs), symvec(ys)))

function _sumout(vars::Vector{Symbol}, expr::IDExpression)
    isempty(vars) ? expr : IDSum(vars, expr)
end

function _product_expr(factors::Vector{IDExpression})
    length(factors) == 1 ? only(factors) : IDProduct(factors)
end

function _remove_incoming(g::ADMG, nodes)
    ns = Set(symvec(nodes))
    ADMG(
        copy(g.vertices),
        [e for e in g.directed_edges if !(e[2] in ns)],
        copy(g.bidirected_edges),
        Dict(k => copy(v) for (k, v) in g.multivariate_variables),
        copy(g.fixed),
    )
end

function _kernel_product(law::IDExpression, vars::Vector{Symbol}, order::Vector{Symbol})
    ordered = [v for v in order if v in vars]
    factors = IDExpression[]
    for v in ordered
        previous = order[1:findfirst(==(v), order)-1]
        push!(factors, IDKernel([v], previous, law))
    end
    _product_expr(factors)
end

function _id_recursive(y::Vector{Symbol}, x::Vector{Symbol},
                       law::IDExpression, g::ADMG)
    V = g.vertices
    Y = _ordered_intersect(y, V)
    X = _ordered_intersect(x, V)

    if isempty(X)
        return _sumout(setdiff(V, Y), law)
    end

    an_y = ancestors(g, Y)
    if !same_set(an_y, V)
        marginalized = _sumout(setdiff(V, an_y), law)
        return _id_recursive(Y, _ordered_intersect(X, an_y),
                             marginalized, subgraph(g, an_y))
    end

    g_xbar = _remove_incoming(g, X)
    an_y_xbar = ancestors(g_xbar, Y)
    W = setdiff(setdiff(V, X), an_y_xbar)
    if !isempty(W)
        return _id_recursive(Y, _ordered_union(X, W), law, g)
    end

    gx = subgraph(g, setdiff(V, X))
    districts_x = all_districts(gx)
    if length(districts_x) > 1
        factors = IDExpression[]
        for S in districts_x
            factor = _id_recursive(S, setdiff(V, S), law, g)
            factor isa _IDFailure && return factor
            push!(factors, factor)
        end
        return _sumout(setdiff(V, _ordered_union(Y, X)), _product_expr(factors))
    end

    S = only(districts_x)
    districts_g = all_districts(g)
    if length(districts_g) == 1
        return _IDFailure((S=S, F=copy(V), treatment=X, outcome=Y))
    end

    order = top_order(g)
    if any(D -> same_set(D, S), districts_g)
        return _sumout(setdiff(S, Y), _kernel_product(law, S, order))
    end

    super_idx = findfirst(D -> all(v -> v in D, S), districts_g)
    superdistrict = super_idx === nothing ? nothing : districts_g[super_idx]
    superdistrict === nothing &&
        return _IDFailure((S=S, F=copy(V), treatment=X, outcome=Y))

    newlaw = _kernel_product(law, superdistrict, order)
    return _id_recursive(Y, _ordered_intersect(X, superdistrict),
                         newlaw, subgraph(g, superdistrict))
end

"""
    ID_algorithm(g::ADMG, treatment, outcome)

Run the Pearl-Shpitser ID algorithm for the ADMG total-effect query
`P(outcome | do(treatment))`.

The return value is an `ADMGIDResult`. If `identified == true`,
`expression` contains a symbolic identifying functional. If identification
fails, `hedge` records the two district sets involved in the failure witness.
"""
function ID_algorithm(g::ADMG, treatment, outcome)
    A = treatment isa AbstractVector ? symvec(treatment) : [sym(treatment)]
    Y = outcome isa AbstractVector ? symvec(outcome) : [sym(outcome)]
    result = _id_recursive(Y, A, IDJoint(copy(g.vertices)), g)
    result isa _IDFailure &&
        return ADMGIDResult(false, A, Y, nothing, result.hedge)
    ADMGIDResult(true, A, Y, result, nothing)
end

ID_algorithm(g::ADMG; treatment, outcome) = ID_algorithm(g, treatment, outcome)

const id_algorithm = ID_algorithm

is_id_identified(g::ADMG, treatment, outcome) =
    ID_algorithm(g, treatment, outcome).identified

summarize_ID(id::ADMGIDResult) =
    (identified=id.identified, treatment=id.treatment, outcome=id.outcome,
     expression=id.expression, hedge=id.hedge)

function Base.show(io::IO, expr::IDJoint)
    print(io, "P(", join(expr.vars, ", "), ")")
end

function Base.show(io::IO, expr::IDKernel)
    print(io, "P(", join(expr.vars, ", "))
    isempty(expr.given) || print(io, " | ", join(expr.given, ", "))
    print(io, ")")
end

function Base.show(io::IO, expr::IDProduct)
    print(io, join(string.(expr.factors), " * "))
end

function Base.show(io::IO, expr::IDSum)
    print(io, "sum_{", join(expr.vars, ", "), "} ", expr.expr)
end

# ── Public fixing/nested helpers ──────────────────────────────────────────────

"""
    fixing_sequence(g, nodes)

Try to fix `nodes` by repeatedly applying the ADMG fixing criterion.
Returns `(fixable, fixing_order, graph)`.
"""
function fixing_sequence(g::ADMG, nodes)
    curr = g
    remaining = nodes isa AbstractVector ? symvec(nodes) : [sym(nodes)]
    order = Symbol[]
    while !isempty(remaining)
        progress = false
        for v in copy(remaining)
            if _is_fixable(curr, v)
                curr = fixed_graph(curr, [v])
                filter!(!=(v), remaining)
                push!(order, v)
                progress = true
                break
            end
        end
        progress || return (fixable=false, fixing_order=order, graph=curr)
    end
    (fixable=true, fixing_order=order, graph=curr)
end

"""
    nested_fixability(g, treatment, outcome)

Run the package's One-Line-ID style nested-fixability check for a total-effect
query. The result includes the SWIG ancestor set `ystar`, the induced graph,
the nested topological order, and per-district reachable-closure diagnostics.
"""
function nested_fixability(g::ADMG, treatment, outcome)
    Aname = sym(treatment)
    Yname = sym(outcome)

    swig = fixed_graph(g, [Aname])
    anc_Y = ancestors(swig, [Yname])
    ystar = [v for v in anc_Y if !get(swig.fixed, v, false)]
    Gystar = subgraph(g, ystar)

    district_results = []
    for D in all_districts(Gystar)
        rc = reachable_closure(g, D)
        push!(district_results, (district=D, reachable_closure=rc.reachable_closure,
                                 fixing_order=rc.fixing_order))
        same_set(rc.reachable_closure, D) ||
            return (identified=false, strategy=:not_identified,
                    treatment=Aname, outcome=Yname, ystar=ystar,
                    Gystar=Gystar, n_order=Symbol[],
                    district_results=district_results)
    end

    n_order = nested_toporder(g, Aname, ystar)
    return (identified=true, strategy=:nested_fixable,
            treatment=Aname, outcome=Yname, ystar=ystar,
            Gystar=Gystar, n_order=n_order,
            district_results=district_results)
end

"""
    is_nested_fixable(g, treatment, outcome)

Return `true` when `nested_fixability(g, treatment, outcome)` identifies the
effect by the package's nested/fixing check.
"""
is_nested_fixable(g::ADMG, treatment, outcome) =
    nested_fixability(g, treatment, outcome).identified

# ── Identification routing ────────────────────────────────────────────────────
#
# Returns a named tuple with at minimum:
#   strategy  — one of :a_fixable, :p_fixable, :nested_fixable, :id_algorithm,
#               :not_identified
#   treatment, outcome  — Symbols
#
# For :nested_fixable additionally:
#   ystar, Gystar, n_order

function identify(g::ADMG, treatment, outcome)
    Aname = sym(treatment)
    Yname = sym(outcome)

    # a-fixable (backdoor): desc(A) ∩ dist(A) = {A}
    if is_fix(g, Aname)
        return (strategy=:a_fixable, treatment=Aname, outcome=Yname)
    end

    # p-fixable (front-door / NPS): dist(A) ∩ ch(A) = {}
    if is_p_fix(g, Aname)
        return (strategy=:p_fixable, treatment=Aname, outcome=Yname)
    end

    nested = nested_fixability(g, Aname, Yname)
    if nested.identified
        return (strategy    = :nested_fixable,
                treatment   = Aname,
                outcome     = Yname,
                ystar       = nested.ystar,
                Gystar      = nested.Gystar,
                n_order     = nested.n_order,
                district_results = nested.district_results)
    end

    id = ID_algorithm(g, Aname, Yname)
    if id.identified
        return (strategy=:id_algorithm, treatment=Aname, outcome=Yname,
                id_expression=id.expression, id_result=id)
    end

    return (strategy=:not_identified, treatment=Aname, outcome=Yname,
            hedge=id.hedge)
end
