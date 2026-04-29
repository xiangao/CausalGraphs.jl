# ── Identification routing ────────────────────────────────────────────────────
#
# Returns a named tuple with at minimum:
#   strategy  — one of :a_fixable, :p_fixable, :nested_fixable, :not_identified
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

    # nested-fixable: One-Line ID
    # Build SWIG: fix treatment vertices
    swig = fixed_graph(g, [Aname])
    # Y* = non-fixed ancestors of outcome in the SWIG
    anc_Y = ancestors(swig, [Yname])
    ystar = [v for v in anc_Y if !get(swig.fixed, v, false)]
    Gystar = subgraph(g, ystar)

    # For each district D in G_Y*, check V \ D is fully fixable in g
    for D in all_districts(Gystar)
        rc = reachable_closure(g, D)
        sort(rc.reachable_closure) == sort(D) ||
            return (strategy=:not_identified, treatment=Aname, outcome=Yname)
    end

    n_order = nested_toporder(g, Aname, ystar)
    return (strategy    = :nested_fixable,
            treatment   = Aname,
            outcome     = Yname,
            ystar       = ystar,
            Gystar      = Gystar,
            n_order     = n_order)
end
