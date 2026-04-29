const Edge = Tuple{Symbol,Symbol}

struct ADMG
    vertices::Vector{Symbol}
    directed_edges::Vector{Edge}
    bidirected_edges::Vector{Edge}
    multivariate_variables::Dict{Symbol,Vector{Symbol}}
    fixed::Dict{Symbol,Bool}
end

# ── Symbol helpers ────────────────────────────────────────────────────────────
sym(x::Symbol) = x
sym(x::AbstractString) = Symbol(x)
symvec(xs) = [sym(x) for x in xs]
edge(x) = (sym(x[1]), sym(x[2]))

function make_graph(; vertices, di_edges=Edge[], bi_edges=Edge[],
                    multivariate_variables=Dict{Symbol,Vector{Symbol}}())
    vs = symvec(vertices)
    mv = Dict{Symbol,Vector{Symbol}}(sym(k) => symvec(v) for (k, v) in pairs(multivariate_variables))
    ADMG(vs, [edge(e) for e in di_edges], [edge(e) for e in bi_edges],
         mv, Dict(v => false for v in vs))
end

function with_graph(graph::Union{ADMG,Nothing}, vertices, di_edges, bi_edges, mv)
    graph === nothing || return graph
    make_graph(vertices=vertices, di_edges=di_edges, bi_edges=bi_edges, multivariate_variables=mv)
end

function insert_vector(vec::Vector{Symbol}, pos::Int, value::Vector{Symbol})
    vcat(vec[1:pos-1], value, vec[pos+1:end])
end

function replace_vector(vec, mv::Dict{Symbol,Vector{Symbol}})
    out = symvec(vec)
    for (var, replacements) in mv
        idx = findfirst(==(var), out)
        idx === nothing || (out = insert_vector(out, idx, replacements))
    end
    out
end

has_biedge(g::ADMG, a::Symbol, b::Symbol) =
    any(e -> (e == (a,b) || e == (b,a)), g.bidirected_edges)
has_diedge(g::ADMG, a::Symbol, b::Symbol) =
    any(e -> e == (a,b), g.directed_edges)

# ── Adjacency & topological sort ─────────────────────────────────────────────
function adjacency_matrix(g::ADMG)
    n = length(g.vertices)
    idx = Dict(v => i for (i,v) in enumerate(g.vertices))
    mat = zeros(Int, n, n)
    for (f,t) in g.directed_edges; mat[idx[t], idx[f]] = 1; end
    mat
end

function _treatment_queue!(queue, vertices, treatment)
    treatment === nothing && return queue
    ti = findfirst(==(sym(treatment)), vertices)
    ti === nothing && return queue
    if ti in queue; filter!(!=(ti), queue); push!(queue, ti); end
    queue
end

function top_order(g::ADMG; treatment=nothing)
    mat = adjacency_matrix(g)
    n = size(mat, 1)
    indeg = vec(sum(mat, dims=2))
    queue = findall(==(0), indeg)
    _treatment_queue!(queue, g.vertices, treatment)
    order = Int[]
    while !isempty(queue)
        node = popfirst!(queue)
        push!(order, node)
        for child in findall(!=(0), mat[:, node])
            indeg[child] -= 1
            if indeg[child] == 0
                push!(queue, child)
                _treatment_queue!(queue, g.vertices, treatment)
            end
        end
    end
    length(order) == n || error("Graph contains a cycle.")
    g.vertices[order]
end

# ── Neighbourhood traversal ───────────────────────────────────────────────────
function parents(g::ADMG, nodes)
    ns = nodes isa AbstractVector ? symvec(nodes) : [sym(nodes)]
    mat = adjacency_matrix(g)
    idx = Dict(v => i for (i,v) in enumerate(g.vertices))
    out = Symbol[]
    for node in ns; append!(out, g.vertices[findall(!=(0), mat[idx[node], :])]); end
    unique(out)
end

function children(g::ADMG, nodes)
    ns = nodes isa AbstractVector ? symvec(nodes) : [sym(nodes)]
    mat = adjacency_matrix(g)
    idx = Dict(v => i for (i,v) in enumerate(g.vertices))
    out = Symbol[]
    for node in ns; append!(out, g.vertices[findall(!=(0), mat[:, idx[node]])]); end
    unique(out)
end

function descendants(g::ADMG, nodes)
    out = nodes isa AbstractVector ? symvec(nodes) : [sym(nodes)]
    frontier = copy(out)
    while !isempty(frontier)
        kids = setdiff(children(g, frontier), out)
        append!(out, kids); frontier = kids
    end
    unique(out)
end

function ancestors(g::ADMG, nodes)
    out = nodes isa AbstractVector ? symvec(nodes) : [sym(nodes)]
    frontier = copy(out)
    while !isempty(frontier)
        pas = setdiff(parents(g, frontier), out)
        append!(out, pas); frontier = pas
    end
    unique(out)
end

function district(g::ADMG, node)
    connected = [sym(node)]
    frontier = copy(connected)
    while !isempty(frontier)
        new_nodes = Symbol[]
        for v in frontier, (a,b) in g.bidirected_edges
            if v == a && !(b in connected); push!(new_nodes, b)
            elseif v == b && !(a in connected); push!(new_nodes, a); end
        end
        append!(connected, unique(new_nodes)); frontier = unique(new_nodes)
    end
    unique(connected)
end

function all_districts(g::ADMG)
    seen = Set{Symbol}()
    result = Vector{Vector{Symbol}}()
    for v in g.vertices
        v in seen && continue
        d = district(g, v)
        push!(result, d)
        union!(seen, d)
    end
    result
end

function markov_blanket(g::ADMG, node)
    n = sym(node)
    d = district(g, n)
    setdiff(unique(vcat(d, parents(g, d))), [n])
end

function markov_pillow(g::ADMG, node; treatment=nothing)
    n = sym(node)
    order = top_order(g; treatment=treatment)
    mb = markov_blanket(g, n)
    upto = findfirst(==(n), order)
    [v for v in mb if v in order[1:upto]]
end

function markov_pillow_with_order(g::ADMG, node::Symbol, order::Vector{Symbol})
    mb = markov_blanket(g, node)
    upto = findfirst(==(node), order)
    upto === nothing && return Symbol[]
    [v for v in mb if v in order[1:upto]]
end

# ── Subgraph ──────────────────────────────────────────────────────────────────
function subgraph(g::ADMG, vertices)
    vs = symvec(vertices)
    vs_set = Set(vs)
    di = [e for e in g.directed_edges  if e[1] in vs_set && e[2] in vs_set]
    bi = [e for e in g.bidirected_edges if e[1] in vs_set && e[2] in vs_set]
    mv = Dict(k => v for (k,v) in g.multivariate_variables if k in vs_set)
    fixed = Dict(v => get(g.fixed, v, false) for v in vs)
    ADMG(vs, di, bi, mv, fixed)
end

# ── Graph modification ────────────────────────────────────────────────────────
function fixed_graph(g::ADMG, nodes)
    out = ADMG(copy(g.vertices), copy(g.directed_edges), copy(g.bidirected_edges),
               Dict(k => copy(v) for (k,v) in g.multivariate_variables),
               copy(g.fixed))
    for n0 in (nodes isa AbstractVector ? symvec(nodes) : [sym(nodes)])
        fixed = copy(out.fixed); fixed[n0] = true
        ps   = parents(out, n0)
        di   = [e for e in out.directed_edges if !(e in [(p,n0) for p in ps])]
        bi   = [e for e in out.bidirected_edges if !(n0 in e)]
        out = ADMG(copy(out.vertices), di, bi,
                   Dict(k => copy(v) for (k,v) in out.multivariate_variables), fixed)
    end
    out
end

# ── Identification primitives ─────────────────────────────────────────────────
function is_fix(g::ADMG, treatment)
    a = sym(treatment)
    length(intersect(descendants(g, [a]), district(g, a))) == 1
end

function is_p_fix(g::ADMG, treatment)
    a = sym(treatment)
    isempty(intersect(district(g, a), children(g, a)))
end

function reachable_closure(g::ADMG, nodes)
    keep = nodes isa AbstractVector ? symvec(nodes) : [sym(nodes)]
    remaining = setdiff(g.vertices, keep)
    order = Symbol[]
    curr = g
    fixed_any = true
    while !isempty(remaining) && fixed_any
        fixed_any = false
        for v in copy(remaining)
            if length(intersect(descendants(curr, v), district(curr, v))) == 1
                curr = fixed_graph(curr, v)
                filter!(!=(v), remaining)
                push!(order, v)
                fixed_any = true
                break
            end
        end
    end
    (reachable_closure=setdiff(g.vertices, order), fixing_order=order, graph=curr)
end

function _is_fixable(g::ADMG, v::Symbol)
    length(intersect(descendants(g, [v]), district(g, v))) == 1
end

function _set_fixable(g::ADMG, vertices::Vector{Symbol})
    curr = g
    remaining = copy(vertices)
    while !isempty(remaining)
        progress = false
        for v in remaining
            if _is_fixable(curr, v)
                curr = fixed_graph(curr, [v])
                filter!(!=(v), remaining)
                progress = true
                break
            end
        end
        progress || return (fixable=false, graph=curr)
    end
    (fixable=true, graph=curr)
end

# ── NPS helpers ────────────────────────────────────────────────────────────────
function is_np_saturated(g::ADMG)
    order = top_order(g)
    for i in 1:length(g.vertices)-1, j in i+1:length(g.vertices)
        vi, vj = g.vertices[i], g.vertices[j]
        v1, v2 = findfirst(==(vi), order) > findfirst(==(vj), order) ? (vj, vi) : (vi, vj)
        rc_v2   = reachable_closure(g, v2).graph
        rc_pair = reachable_closure(g, [vi, vj]).graph
        if !(v1 in parents(g, district(rc_v2, v2))) &&
           (length(all_districts(rc_pair)) - length(reachable_closure(g, [vi, vj]).fixing_order) > 1)
            return false
        end
    end
    true
end

function cml(g::ADMG, treatment)
    a = sym(treatment)
    order = top_order(g; treatment=a)
    apos = findfirst(==(a), order)
    C = order[1:apos-1]
    L = intersect(district(g, a), order[apos:end])
    M = setdiff(g.vertices, vcat(C, L))
    (C=C, M=M, L=L)
end

# ── Nested-fixable topology ───────────────────────────────────────────────────
function nested_toporder(g::ADMG, Aname::Symbol, ystar::Vector{Symbol})
    dist_T = district(g, Aname)
    focal = intersect(ystar, dist_T)
    isempty(focal) && (focal = [Aname])
    desc_focal = descendants(g, focal)
    non_desc   = setdiff(g.vertices, desc_focal)
    vcat(top_order(subgraph(g, non_desc)), top_order(subgraph(g, desc_focal)))
end

# ── MB-shielded ────────────────────────────────────────────────────────────────
function mb_shielded(g::ADMG)
    vs = g.vertices
    n  = length(vs)
    for i in 1:n-1, j in i+1:n
        vi, vj = vs[i], vs[j]
        adjacent = has_biedge(g, vi, vj) || has_diedge(g, vi, vj) || has_diedge(g, vj, vi)
        if !adjacent && (vi in markov_blanket(g, vj) || vj in markov_blanket(g, vi))
            return false
        end
    end
    true
end

# ── Scalar helpers ────────────────────────────────────────────────────────────
rerank(target, reference) = sort(symvec(target), by=x -> findfirst(==(x), reference))
same_set(a, b) = sort(symvec(a)) == sort(symvec(b))

col(df::DataFrame, name) = Float64.(df[!, sym(name)])
nrows(df::DataFrame) = nrow(df)
is_binary(x) = all(v -> v == 0 || v == 1, skipmissing(x))
clip(x, lo, hi) = min.(max.(x, lo), hi)
safe_logit(p) = log.(clip(p, 1e-8, 1-1e-8) ./ (1 .- clip(p, 1e-8, 1-1e-8)))
expit(x) = 1.0 ./ (1.0 .+ exp.(-x))
probability_a(pA1, a) = a == 1 ? pA1 : 1 .- pA1
normal_pdf(x, μ, σ) = exp(-0.5*((x-μ)/max(Float64(σ),1e-8))^2) / (sqrt(2π)*max(Float64(σ),1e-8))

product_columns(dict, vars, n) =
    isempty(vars) ? ones(n) : reduce(.*, [dict[v] for v in vars])

function set_column(df::DataFrame, name::Symbol, value)
    out = copy(df); out[!, name] = fill(value, nrow(out)); out
end

function ci_from_eif(estimate, eif, n)
    se = sqrt(mean(eif .^ 2) / n)
    (estimate - 1.96*se, estimate + 1.96*se)
end

# ── Nuisance math helpers ─────────────────────────────────────────────────────
function design_matrix(df::DataFrame, predictors::Vector{Symbol})
    n = nrow(df)
    X = ones(Float64, n, length(predictors)+1)
    for (j,p) in enumerate(predictors); X[:,j+1] = col(df, p); end
    X
end

function solve_ridge(X, y; weights=nothing, λ=1e-8)
    if weights === nothing
        (X'X + λ*I(size(X,2))) \ X'y
    else
        w = Float64.(weights)
        ((X' .* w') * X + λ*I(size(X,2))) \ (X' * (y .* w))
    end
end

fit_lm(df::DataFrame, y, predictors::Vector{Symbol}; weights=nothing) =
    solve_ridge(design_matrix(df, predictors), Float64.(y); weights=weights)
predict_lm(β, df::DataFrame, predictors::Vector{Symbol}) =
    design_matrix(df, predictors) * β

function fit_logistic(df::DataFrame, y, predictors::Vector{Symbol};
                       maxiter=100, tol=1e-8, sample_weights=nothing)
    X  = design_matrix(df, predictors)
    yy = Float64.(y)
    sw = sample_weights === nothing ? ones(size(X,1)) : Float64.(sample_weights)
    β  = zeros(size(X,2))
    for _ in 1:maxiter
        η   = X * β
        μ   = clip(expit(η), 1e-8, 1-1e-8)
        W   = μ .* (1 .- μ) .* sw
        z   = η + (yy - μ) ./ max.(μ .* (1 .- μ), 1e-10)
        βnew = solve_ridge(X, z; weights=W, λ=1e-7)
        norm(βnew - β) < tol && return βnew
        β = βnew
    end
    β
end

predict_logistic(β, df::DataFrame, predictors::Vector{Symbol}) =
    clip(expit(design_matrix(df, predictors) * β), 1e-8, 1-1e-8)

function scalar_logistic_fluctuation(y, offset, H; maxiter=100, tol=1e-9)
    yy, off, hh = Float64.(y), Float64.(offset), Float64.(H)
    ε = 0.0
    for _ in 1:maxiter
        μ = clip(expit(off .+ ε .* hh), 1e-8, 1-1e-8)
        score = sum(hh .* (yy .- μ))
        info  = sum((hh .^ 2) .* μ .* (1 .- μ))
        info <= eps() && return ε
        step = score / info; ε += step
        abs(step) < tol && return ε
    end
    ε
end

weighted_mean_residual(y, offset, w) =
    sum(Float64.(w) .* (Float64.(y) .- Float64.(offset))) / max(sum(Float64.(w)), eps())
