const MEdge = Tuple{String,String}

struct MDAG
    vertices::Vector{String}
    fixed::Dict{String,Bool}
    di_edges::Vector{MEdge}
    obs_variables::Vector{String}
    missing_variables::Vector{String}
    missing_indicators::Vector{String}
end

mutable struct MissingTree
    root::String
    leafs::Dict{String,MissingTree}
    di_edges::Vector{MEdge}
    pruned_nodes::Dict{String,Bool}
    pruned_edges::Dict{MEdge,Bool}
end

struct IDLaw
    graph::MDAG
    id_status::Bool
    id_status_collection::Dict{String,Bool}
    tree_collection::Dict{String,MissingTree}
    selection_collection::Dict{String,Vector{String}}
    selectionx_collection::Dict{String,Vector{String}}
    selectionr_collection::Dict{String,Vector{String}}
    preselection_collection::Dict{String,Vector{String}}
    is_fullID::Bool
    id_targetlaw::Union{Nothing,IDLaw}
end

struct IDResult
    target::IDLaw
    full::Union{Nothing,IDLaw}
end

struct PropensityEstimate
    fit::Any
    pred::Vector{Union{Missing,Float64}}
    parents::Vector{String}
    pred_rows::Vector{Bool}
    fit_rows::Vector{Bool}
    selection_set::Vector{String}
    free_R::Vector{String}
end

struct CliqueResult
    clique::Vector{String}
    magic_weight::Vector{Float64}
    fit_rows::Vector{Bool}
end

struct MEstimationResult
    estimates::Vector{Float64}
    vcov::Matrix{Float64}
    se::Vector{Float64}
    converged::Bool
    iterations::Int
    estimating_equations::Vector{Float64}
end

_vec(x) = String[string(v) for v in x]
_edge(e) = (string(e[1]), string(e[2]))
_edges(edges) = MEdge[_edge(e) for e in edges]

function _unique(xs)
    out = String[]
    seen = Set{String}()
    for x in xs
        sx = string(x)
        if !(sx in seen)
            push!(out, sx)
            push!(seen, sx)
        end
    end
    out
end

_intersect(a, b) = [x for x in a if x in Set(b)]
_setdiff(a, b) = [x for x in a if !(x in Set(b))]
_union(a, b) = _unique(vcat(a, b))

function make_mdag(; obs_variables=String[], missing_variables=String[], missing_indicators=String[], di_edges)
    obs = _vec(obs_variables)
    missx = _vec(missing_variables)
    missr = _vec(missing_indicators)
    vertices = vcat(obs, missx, missr)
    isempty(vertices) && error("At least one of obs_variables or missing_variables must be non-empty.")
    length(missx) == length(missr) || error("The length of missing_variables and missing_indicators must be the same.")
    fixed = Dict(v => false for v in vertices)
    MDAG(vertices, fixed, _edges(di_edges), obs, missx, missr)
end

function make_mtree(root; leafs=Dict{String,MissingTree}(), di_edges=MEdge[],
                    pruned_nodes=nothing, pruned_edges=nothing)
    r = string(root)
    leaves = Dict(string(k) => deepcopy(v) for (k, v) in leafs)
    edges = MEdge[_edge(e) for e in di_edges]
    pn = pruned_nodes === nothing ? Dict(k => false for k in keys(leaves)) :
         Dict(string(k) => Bool(v) for (k, v) in pruned_nodes)
    pe = if pruned_edges === nothing
        Dict(e => get(pn, e[2], false) for e in edges)
    else
        Dict(_edge(k) => Bool(v) for (k, v) in pruned_edges)
    end
    MissingTree(r, leaves, edges, pn, pe)
end

function append_leaf(Tk::MissingTree, Ti::MissingTree; pruned=false)
    out = deepcopy(Tk)
    child = deepcopy(Ti)
    if !(child.root in keys(out.leafs))
        out.pruned_nodes[child.root] = pruned
        e = (out.root, child.root)
        push!(out.di_edges, e)
        out.pruned_edges[e] = pruned
    end
    out.leafs[child.root] = child
    out
end

function prune_tree(tree::MissingTree, nodes)
    out = deepcopy(tree)
    to_prune = Set(_vec(nodes isa AbstractString ? [nodes] : nodes))
    for n in to_prune
        out.pruned_nodes[n] = true
    end
    for e in out.di_edges
        out.pruned_edges[e] = get(out.pruned_nodes, e[2], false)
    end
    out
end

function get_all_nodes(tree::MissingTree)
    nodes = [tree.root]
    for child in values(tree.leafs)
        append!(nodes, get_all_nodes(child))
    end
    _unique(nodes)
end

function get_all_edges(tree::MissingTree)
    edges = copy(tree.di_edges)
    pruned = copy(tree.pruned_edges)
    for child in values(tree.leafs)
        ce, cp = get_all_edges(child)
        append!(edges, ce)
        merge!(pruned, cp)
    end
    edges, pruned
end

function are_trees_identical(a::MissingTree, b::MissingTree)
    a.root == b.root || return false
    a.di_edges == b.di_edges || return false
    a.pruned_nodes == b.pruned_nodes || return false
    a.pruned_edges == b.pruned_edges || return false
    sort(collect(keys(a.leafs))) == sort(collect(keys(b.leafs))) || return false
    all(are_trees_identical(a.leafs[k], b.leafs[k]) for k in keys(a.leafs))
end

function adj_matrix(graph::MDAG)
    n = length(graph.vertices)
    idx = Dict(v => i for (i, v) in pairs(graph.vertices))
    mat = zeros(Int, n, n)
    for (from, to) in graph.di_edges
        mat[idx[to], idx[from]] = 1
    end
    mat
end

function _treatment_queue(vertices, queue, treatment)
    treatment === nothing && return queue
    ti = findfirst(==(string(treatment)), vertices)
    if ti !== nothing && ti in queue
        queue = [q for q in queue if q != ti]
        push!(queue, ti)
    end
    queue
end

function _top_order_mat(mat, vertices; treatment=nothing)
    n = size(mat, 1)
    indeg = vec(sum(mat, dims=2))
    queue = findall(==(0), indeg)
    queue = _treatment_queue(vertices, queue, treatment)
    order = Int[]
    while !isempty(queue)
        node = popfirst!(queue)
        push!(order, node)
        for child in findall(!=(0), mat[:, node])
            indeg[child] -= 1
            if indeg[child] == 0
                push!(queue, child)
                queue = _treatment_queue(vertices, queue, treatment)
            end
        end
    end
    length(order) == n || error("The graph contains a cycle. Topological ordering is only available for DAGs.")
    vertices[order]
end

top_order(graph::MDAG; treatment=nothing) = _top_order_mat(adj_matrix(graph), graph.vertices; treatment)
top_order_R(graph::MDAG; treatment=nothing) = [v for v in top_order(graph; treatment) if v in graph.missing_indicators]

function parents(graph::MDAG, nodes)
    ns = _vec(nodes isa AbstractString ? [nodes] : nodes)
    mat = adj_matrix(graph)
    idx = Dict(v => i for (i, v) in pairs(graph.vertices))
    out = String[]
    for n in ns
        append!(out, graph.vertices[findall(!=(0), mat[idx[n], :])])
    end
    _unique(out)
end

function children(graph::MDAG, nodes)
    ns = _vec(nodes isa AbstractString ? [nodes] : nodes)
    mat = adj_matrix(graph)
    idx = Dict(v => i for (i, v) in pairs(graph.vertices))
    out = String[]
    for n in ns
        append!(out, graph.vertices[findall(!=(0), mat[:, idx[n]])])
    end
    _unique(out)
end

function children(tree::MissingTree, node)
    string(node) == tree.root || error("When the input graph is a tree, node must be the root node.")
    [e[2] for e in tree.di_edges if !get(tree.pruned_nodes, e[2], false)]
end

function descendants(graph::MDAG, nodes)
    out = _vec(nodes isa AbstractString ? [nodes] : nodes)
    frontier = copy(out)
    while !isempty(frontier)
        ch = children(graph, frontier)
        new = _setdiff(ch, out)
        append!(out, new)
        frontier = new
    end
    out
end

function ancestors(graph::MDAG, nodes)
    out = _vec(nodes isa AbstractString ? [nodes] : nodes)
    frontier = copy(out)
    while !isempty(frontier)
        pa = parents(graph, frontier)
        new = _setdiff(pa, out)
        append!(out, new)
        frontier = new
    end
    out
end

nondescendants(graph::MDAG, nodes) = _setdiff(graph.vertices, descendants(graph, nodes))

function gfix(graph::MDAG, nodes)
    out = deepcopy(graph)
    ns = _vec(nodes isa AbstractString ? [nodes] : nodes)
    for n in ns
        out.fixed[n] = true
        ps = parents(out, n)
        remove = Set((p, n) for p in ps)
        filter!(e -> !(e in remove), out.di_edges)
    end
    out
end

function r_indicator_for_x(graph::MDAG, xs)
    set = _vec(xs isa AbstractString ? [xs] : xs)
    all(x -> x in graph.missing_variables, set) || error("Input must contain only missing variables.")
    [graph.missing_indicators[findfirst(==(x), graph.missing_variables)] for x in set]
end

function x_variable_for_r(graph::MDAG, rs)
    set = _vec(rs isa AbstractString ? [rs] : rs)
    all(r -> r in graph.missing_indicators, set) || error("Input must contain only missing indicators.")
    [graph.missing_variables[findfirst(==(r), graph.missing_indicators)] for r in set]
end

function _findRX(R, X, set)
    s = _vec(set isa AbstractString ? [set] : set)
    if all(x -> x in R, s)
        return [X[findfirst(==(x), R)] for x in s]
    elseif all(x -> x in X, s)
        return [R[findfirst(==(x), X)] for x in s]
    end
    error("The input set needs to be all missing indicators or all missing variables.")
end

function selectionx_set(graph::MDAG, nodes)
    ns = _vec(nodes isa AbstractString ? [nodes] : nodes)
    out = String[]
    for n in ns
        ps = parents(graph, n)
        append!(out, graph.missing_indicators[findall(x -> x in ps, graph.missing_variables)])
    end
    _unique(out)
end

function problematic_set(graph::MDAG, nodes)
    _intersect(descendants(graph, nodes), selectionx_set(graph, nodes))
end

function colluder(graph::MDAG, Rk; Rj=nothing)
    Rj = Rj === nothing ? Rk : Rj
    Rks = string(Rk)
    Rjs = string(Rj)
    Rks in graph.missing_indicators || error("Rk must be a missing indicator.")
    tau = top_order_R(graph)
    Xj = only(_findRX(graph.missing_indicators, graph.missing_variables, [Rjs]))
    out = String[]
    for d in descendants(graph, Rks)
        if Xj in parents(graph, d)
            push!(out, d)
        end
    end
    out[sortperm([findfirst(==(x), tau) for x in out])]
end

function colluder_chain(graph::MDAG, node)
    tau = top_order_R(graph)
    c = colluder(graph, string(node))
    out = copy(c)
    for rj in _setdiff(c, [string(node)])
        append!(out, colluder_chain(graph, rj))
    end
    out = _unique(out)
    out[sortperm([findfirst(==(x), tau) for x in out])]
end

function _has_directed_edge(graph::MDAG, from, to)
    (string(from), string(to)) in graph.di_edges
end

function _neighbors(graph::MDAG, node)
    out = String[]
    for (a, b) in graph.di_edges
        a == node && push!(out, b)
        b == node && push!(out, a)
    end
    _unique(out)
end

function _all_simple_paths(graph::MDAG, A, B)
    A == B && return Vector{Vector{String}}()
    paths = Vector{Vector{String}}()
    function dfs(node, target, visited, path)
        if node == target
            push!(paths, copy(path))
            return
        end
        for nb in _neighbors(graph, node)
            if !(nb in visited)
                push!(visited, nb)
                push!(path, nb)
                dfs(nb, target, visited, path)
                pop!(path)
                delete!(visited, nb)
            end
        end
    end
    dfs(string(A), string(B), Set([string(A)]), [string(A)])
    paths
end

function _is_collider(graph::MDAG, path::Vector{String}, i::Int)
    prev, node, nxt = path[i - 1], path[i], path[i + 1]
    prev_arrow = _has_directed_edge(graph, prev, node) ? "->" : "<-"
    next_arrow = _has_directed_edge(graph, node, nxt) ? "->" : "<-"
    prev_arrow == "->" && next_arrow == "<-"
end

function is_d_separated(graph::MDAG, A, B, C=String[])
    As, Bs = string(A), string(B)
    Cs = Set(_vec(C isa AbstractString ? [C] : C))
    As == Bs && return (d_separated=false, paths=Vector{Vector{String}}())
    paths = _all_simple_paths(graph, As, Bs)
    isempty(paths) && return (d_separated=true, paths=paths)
    for path in paths
        length(path) == 2 && return (d_separated=false, paths=paths)
    end
    blocked = Bool[]
    for path in paths
        path_blocked = false
        for i in 2:(length(path)-1)
            node = path[i]
            if _is_collider(graph, path, i)
                if !any(d -> d in Cs, descendants(graph, node))
                    path_blocked = true
                    break
                end
            else
                if node in Cs
                    path_blocked = true
                    break
                end
            end
        end
        push!(blocked, path_blocked)
    end
    (d_separated=all(blocked), paths=paths)
end

function _id_status(graph, selection_collection, Rk, preselection_set, Tk; vocal=false)
    prob = problematic_set(graph, Rk)
    nd = nondescendants(graph, Rk)
    pa = parents(graph, Rk)
    ch = children(Tk, Rk)
    fixed_graph = gfix(graph, ch)
    check = _setdiff(preselection_set, nd)
    not_d_sep = String[]
    for x in check
        dsep = x in ch ? true : is_d_separated(fixed_graph, Rk, x, pa).d_separated
        dsep || push!(not_d_sep, x)
    end
    if isempty(not_d_sep)
        return (id_status=true, preselection_set=_unique(preselection_set), tree=Tk)
    elseif !isempty(_intersect(not_d_sep, prob))
        vocal && @info "$Rk is not d-separated from problematic nodes; ID failed."
        return (id_status=false, preselection_set=_unique(preselection_set), tree=Tk)
    end

    tree_collection = Dict{String,MissingTree}()
    local_selection = Dict{String,Vector{String}}()
    Ck = _unique(reduce(vcat, [colluder(graph, Rk; Rj=rj) for rj in not_d_sep]; init=String[]))
    outTk = Tk
    for Ri in ch
        Ti = outTk.leafs[Ri]
        Rich = children(Ti, Ri)
        if Ri in Ck
            outTk = prune_tree(outTk, [Ri])
        else
            condition1 = !isempty(_intersect(Rich, collect(keys(tree_collection))))
            Rm_select_Rj = String[]
            for Rm in Rich
                for Rj in not_d_sep
                    if Rm in colluder(graph, Rk; Rj) && Rj in parents(graph, Ri)
                        push!(Rm_select_Rj, Rm)
                    end
                end
            end
            Rm_select_Rj = _unique(Rm_select_Rj)
            if condition1 || !isempty(Rm_select_Rj)
                result = _tree_prune(graph, selection_collection, Rm_select_Rj, Ri, Ti, Rk,
                                     outTk, Ck, tree_collection, local_selection)
                outTk = result.Tk
                tree_collection = result.tree_collection
                local_selection = result.selection_collection
            end
        end
    end

    sx = selectionx_set(graph, Rk)
    new_pre = copy(sx)
    for Ri in children(outTk, Rk)
        append!(new_pre, haskey(local_selection, Ri) ? local_selection[Ri] : selection_collection[Ri])
    end
    _id_status(graph, selection_collection, Rk, _unique(new_pre), outTk; vocal=false)
end

function _tree_prune(graph, selection_collection, Rm_select_Rj, Ri, Ti, Rk, Tk, Ck,
                     tree_collection, local_selection)
    outTi = Ti
    outTk = Tk
    pre = selectionx_set(graph, Ri)
    for Rm in children(outTi, Ri)
        if haskey(tree_collection, Rm)
            outTi = append_leaf(outTi, tree_collection[Rm])
            append!(pre, local_selection[Rm])
        elseif Rm in Rm_select_Rj
            outTi = prune_tree(outTi, [Rm])
        else
            append!(pre, selection_collection[Rm])
        end
    end
    st = _id_status(graph, selection_collection, Ri, _unique(pre), outTi; vocal=false)
    if st.id_status
        outTk = append_leaf(outTk, st.tree)
        tree_collection[Ri] = st.tree
        sr = _intersect(parents(graph, Ri), st.preselection_set)
        local_selection[Ri] = _union(selectionx_set(graph, Ri), sr)
    else
        outTk = append_leaf(outTk, st.tree)
        outTk = prune_tree(outTk, [Ri])
    end
    (Ti=outTi, Tk=outTk, tree_collection=tree_collection, selection_collection=local_selection)
end

function _tree_construction(graph, Rk, Tk, sx, tree_collection, selection_collection, selectionr_collection, id_status_collection)
    tau = reverse(top_order_R(graph))
    chain = colluder_chain(graph, Rk)
    desc = _setdiff(descendants(graph, Rk), [Rk])
    noid = [k for (k, v) in id_status_collection if !v]
    candidates = _setdiff(desc, vcat(chain, noid))
    candidates = candidates[sortperm([findfirst(==(x), tau) for x in candidates])]
    pre = copy(sx)
    outTk = Tk
    for Ri in candidates
        outTk = append_leaf(outTk, tree_collection[Ri])
        append!(pre, selection_collection[Ri])
    end
    pre = _unique(pre)
    st = _id_status(graph, selection_collection, Rk, pre, outTk)
    pre = st.preselection_set
    sr = _intersect(parents(graph, Rk), pre)
    sel = _union(sx, sr)
    if !st.id_status
        @info "Identification for $Rk failed. The target law is NOT identified in the given mDAG."
    end
    (id_status=st.id_status, tree=st.tree, selection_set=sel, selectionr_set=sr, preselection_set=pre)
end

function ID_algorithm(graph::MDAG; fulllaw=false)
    tau = reverse(top_order_R(graph))
    trees = Dict{String,MissingTree}()
    pre = Dict{String,Vector{String}}()
    sel = Dict{String,Vector{String}}()
    selx = Dict{String,Vector{String}}()
    selr = Dict{String,Vector{String}}()
    ids = Dict{String,Bool}()

    for Rk in tau
        sx = selectionx_set(graph, Rk)
        selx[Rk] = sx
        prob = problematic_set(graph, Rk)
        Tk = make_mtree(Rk)
        if isempty(prob)
            trees[Rk] = Tk
            sr = _intersect(parents(graph, Rk), sx)
            selr[Rk] = sr
            sel[Rk] = _union(sx, sr)
            pre[Rk] = sx
            ids[Rk] = true
        elseif Rk in prob
            @info "Missing mechanism for $Rk is self-censored. The target law is NOT identified for it."
            trees[Rk] = Tk
            sr = _intersect(parents(graph, Rk), sx)
            selr[Rk] = sr
            sel[Rk] = _union(sx, sr)
            pre[Rk] = sx
            ids[Rk] = false
        else
            tc = _tree_construction(graph, Rk, Tk, sx, trees, sel, selr, ids)
            trees[Rk] = tc.tree
            selr[Rk] = tc.selectionr_set
            sel[Rk] = tc.selection_set
            pre[Rk] = tc.preselection_set
            ids[Rk] = tc.id_status
        end
    end

    ok = all(values(ids))
    ok && @info "The target law is identified in the given mDAG."
    target = IDLaw(graph, ok, ids, trees, sel, selx, selr, pre, false, nothing)
    full = fulllaw ? ID_flex_eval(graph, target, selr) : nothing
    if fulllaw
        full.id_status ? @info("The full law is identified in the given mDAG.") :
                         @info("The full law is NOT identified in the given mDAG.")
    end
    IDResult(target, full)
end

function ID_flex_eval(graph::MDAG, id_targetlaw::IDLaw, randomR)
    tau = reverse(top_order_R(graph))
    target_ids = id_targetlaw.id_status_collection
    candidates = [k for (k, v) in target_ids if v]
    candidates = candidates[sortperm([findfirst(==(x), tau) for x in candidates])]
    ids = copy(target_ids)
    sel = deepcopy(id_targetlaw.selection_collection)
    pre = deepcopy(id_targetlaw.preselection_collection)
    selr = deepcopy(id_targetlaw.selectionr_collection)
    selx = deepcopy(id_targetlaw.selectionx_collection)
    trees = deepcopy(id_targetlaw.tree_collection)
    is_full = randomR == id_targetlaw.selectionr_collection

    for Rk in candidates
        Tk = trees[Rk]
        ch = children(Tk, Rk)
        ch = ch[sortperm([findfirst(==(x), tau) for x in ch])]
        rR = get(randomR, Rk, String[])
        if !isempty(_intersect(selr[Rk], rR))
            if isempty(ch) || !isempty(_intersect(rR, selx[Rk]))
                ids[Rk] = false
                @info (is_full ? "Full law identification for $Rk failed." :
                                 "Identification for $Rk at specified evaluation failed.")
            else
                local_trees = Dict{String,MissingTree}()
                local_sel = Dict{String,Vector{String}}()
                Ck = _unique(reduce(vcat, [colluder(graph, Rk; Rj=rj) for rj in rR]; init=String[]))
                for Ri in ch
                    Ti = Tk.leafs[Ri]
                    Rich = children(Ti, Ri)
                    if Ri in Ck
                        Tk = prune_tree(Tk, [Ri])
                    else
                        Rm_select_Rj = String[]
                        for Rm in Rich, Rj in rR
                            if Rm in colluder(graph, Rk; Rj) && Rj in parents(graph, Rk)
                                push!(Rm_select_Rj, Rm)
                            end
                        end
                        Rm_select_Rj = _unique(Rm_select_Rj)
                        if !isempty(Rm_select_Rj)
                            tp = _tree_prune(graph, sel, Rm_select_Rj, Ri, Ti, Rk, Tk, Ck, local_trees, local_sel)
                            Tk = tp.Tk
                            local_trees = tp.tree_collection
                            local_sel = tp.selection_collection
                        end
                    end
                end
                sx = selectionx_set(graph, Rk)
                newpre = copy(sx)
                for Ri in children(Tk, Rk)
                    append!(newpre, haskey(local_sel, Ri) ? local_sel[Ri] : sel[Ri])
                end
                st = _id_status(graph, sel, Rk, _unique(newpre), Tk)
                if !st.id_status
                    ids[Rk] = false
                    @info (is_full ? "Full law identification for $Rk failed." :
                                     "Identification for $Rk at specified evaluation failed.")
                end
                pre[Rk] = st.preselection_set
                selr[Rk] = _intersect(parents(graph, Rk), st.preselection_set)
                sel[Rk] = _union(sx, selr[Rk])
                trees[Rk] = st.tree
            end
        end
    end
    ok = id_targetlaw.id_status ? all(values(ids)) : false
    IDLaw(graph, ok, ids, trees, sel, selx, selr, pre, is_full, id_targetlaw)
end

function _law(ID::IDResult, law)
    law == :target || law == "target" ? ID.target :
    law == :full || law == "full" ? (ID.full === nothing ? error("ID result has no full-law component.") : ID.full) :
    error("law must be :target or :full")
end

function summarize_ID(ID::IDResult)
    rows_for(id_obj) = begin
        rows = NamedTuple[]
        vars = collect(keys(id_obj.id_status_collection))
        tau = reverse(top_order_R(id_obj.graph))
        vars = vars[sortperm([findfirst(==(x), tau) for x in vars])]
        for v in vars
            pa = parents(id_obj.graph, v)
            push!(rows, (
                Missing_Indicators=v,
                Propensity="p($v | $(join(pa, ", ")))",
                ID_Status=id_obj.id_status_collection[v],
                Intervened_On=isempty(children(id_obj.tree_collection[v], v)) ? "-" : join(children(id_obj.tree_collection[v], v), ", "),
                Preselection=isempty(id_obj.preselection_collection[v]) ? "-" : join(id_obj.preselection_collection[v], ", "),
                Selection=isempty(id_obj.selection_collection[v]) ? "-" : join(id_obj.selection_collection[v], ", "),
                Selection_r=isempty(id_obj.selectionr_collection[v]) ? "-" : join(id_obj.selectionr_collection[v], ", "),
                Selection_x=isempty(id_obj.selectionx_collection[v]) ? "-" : join(id_obj.selectionx_collection[v], ", "),
            ))
        end
        DataFrame(rows)
    end
    target = rows_for(ID.target)
    full = ID.full === nothing ? nothing : rows_for(ID.full)
    (target=target, full=full)
end

function _col(data::DataFrame, name::String)
    data[!, name]
end

function _all_one_rows(data::DataFrame, vars::Vector{String})
    isempty(vars) && return trues(nrow(data))
    rows = trues(nrow(data))
    for v in vars
        rows .&= [isequal(x, 1) for x in _col(data, v)]
    end
    rows
end

function _formula(response::String, predictors::Vector{String}, formula)
    if formula === nothing
        rhs = isempty(predictors) ? ConstantTerm(1) : sum(term.(Symbol.(predictors)))
        return term(Symbol(response)) ~ rhs
    end
    formula isa Dict && return formula[response]
    formula
end

function _predict_glm(fit, newdata::DataFrame)
    Float64.(predict(fit, newdata))
end

function _fit_associational(graph, data::DataFrame, Rk, preselection_set, selection_set, selectionr_set;
                            formula=nothing, link=LogitLink(), ml_model=nothing)
    pa = parents(graph, Rk)
    pred_rows = _all_one_rows(data, selection_set)
    fit_rows = _all_one_rows(data, preselection_set)
    data_obs = data[fit_rows, vcat([Rk], pa)]
    free_R = _setdiff(_intersect(pa, graph.missing_indicators), selectionr_set)
    pred = Vector{Union{Missing,Float64}}(missing, nrow(data))
    if ml_model !== nothing && !isempty(pa)
        X_fit = data_obs[:, pa]
        y_fit = Float64.(data_obs[!, Rk])
        mach = _mlj_fit_classifier(ml_model, X_fit, y_fit)
        X_pred = data[pred_rows, pa]
        pred[pred_rows] = _mlj_predict_prob1(mach, X_pred)
        return PropensityEstimate(mach, pred, pa, pred_rows, fit_rows, selection_set, free_R)
    end
    fit = glm(_formula(Rk, pa, formula), data_obs, Binomial(), link)
    # StatsModels.predict(fit, newdata) fails for intercept-only models when
    # newdata has no RHS columns; fall back to predict(fit) on training data.
    pred[pred_rows] = if isempty(pa) && formula === nothing
        p = Float64.(predict(fit))
        fill(p[1], sum(pred_rows))
    else
        _predict_glm(fit, data[pred_rows, vcat([Rk], pa)])
    end
    PropensityEstimate(fit, pred, pa, pred_rows, fit_rows, selection_set, free_R)
end

function _fit_weighted_R(graph, data::DataFrame, weights_collection, Tk, Rk, preselection_set, selection_set, selectionr_set;
                         formula=nothing, link=LogitLink())
    pa = parents(graph, Rk)
    pred_rows = _all_one_rows(data, selection_set)
    fit_rows = _all_one_rows(data, preselection_set)
    free_R = _setdiff(_intersect(pa, graph.missing_indicators), selectionr_set)
    weights = ones(Float64, nrow(data))
    for (pi, est) in weights_collection
        p = [ismissing(x) ? 1.0 : Float64(x) for x in est.pred]
        if pi in preselection_set
            weights .*= p
        else
            rcol = Float64.(_col(data, pi))
            weights .*= rcol .* p .+ (1 .- rcol) .* (1 .- p)
        end
    end
    data_obs = data[fit_rows, vcat([Rk], pa)]
    fit = glm(_formula(Rk, pa, formula), data_obs, Binomial(), link; wts=1 ./ weights[fit_rows])
    pred = Vector{Union{Missing,Float64}}(missing, nrow(data))
    pred[pred_rows] = _predict_glm(fit, data[pred_rows, vcat([Rk], pa)])
    PropensityEstimate(fit, pred, pa, pred_rows, fit_rows, selection_set, free_R)
end

function _estimate_causal(graph, data::DataFrame, tree_collection, Rk, preselection_collection,
                          selection_collection, selectionr_collection, propensitys; formula=nothing, link=LogitLink())
    Tk = tree_collection[Rk]
    ch = children(Tk, Rk)
    weights_collection = Dict{String,PropensityEstimate}()
    local_pre = Dict{String,Vector{String}}()
    local_sel = Dict{String,Vector{String}}()
    local_selr = Dict{String,Vector{String}}()
    changed = String[]

    for Ri in ch
        Ti = Tk.leafs[Ri]
        if are_trees_identical(Ti, tree_collection[Ri])
            weights_collection[Ri] = propensitys[Ri]
        else
            push!(changed, Ri)
            Rich = children(Ti, Ri)
            ripre = selectionx_set(graph, Ri)
            for Rj in Rich
                append!(ripre, Rj in changed ? local_sel[Rj] : selection_collection[Rj])
            end
            ripre = _unique(ripre)
            riselr = _intersect(parents(graph, Ri), ripre)
            risel = _union(selectionx_set(graph, Ri), riselr)
            local_pre[Ri], local_sel[Ri], local_selr[Ri] = ripre, risel, riselr
            child_weights = Dict{String,PropensityEstimate}(Rj => weights_collection[Rj] for Rj in Rich)
            weights_collection[Ri] = _fit_weighted_R(graph, data, child_weights, Ti, Ri, ripre, risel, riselr;
                                                     formula, link)
        end
    end
    _fit_weighted_R(graph, data, weights_collection, Tk, Rk, preselection_collection[Rk],
                    selection_collection[Rk], selectionr_collection[Rk]; formula, link)
end

function propensity(graph::MDAG, data::DataFrame, ID::IDResult; law=:target, formula=nothing, link=LogitLink(), vocal=true, model=nothing)
    idlaw = _law(ID, law)
    nodes = [k for (k, v) in idlaw.id_status_collection if v]
    tau = reverse(top_order_R(graph))
    nodes = nodes[sortperm([findfirst(==(x), tau) for x in nodes])]
    out = Dict{String,PropensityEstimate}()
    for Rk in nodes
        vocal && @info "Estimating propensity for $Rk."
        Tk = idlaw.tree_collection[Rk]
        if isempty(children(Tk, Rk))
            out[Rk] = _fit_associational(graph, data, Rk, idlaw.preselection_collection[Rk],
                                         idlaw.selection_collection[Rk], idlaw.selectionr_collection[Rk];
                                         formula, link, ml_model=model)
        else
            out[Rk] = _estimate_causal(graph, data, idlaw.tree_collection, Rk, idlaw.preselection_collection,
                                       idlaw.selection_collection, idlaw.selectionr_collection, out;
                                       formula, link)
        end
    end
    out
end

function find_clique(graph::MDAG, data::DataFrame, propensitys::Dict{String,PropensityEstimate}, est_R, non_ID_nodes=String[]; vocal=true)
    clique = _vec(est_R isa AbstractString ? [est_R] : est_R)
    processed = String[]
    while !isempty(_setdiff(clique, processed))
        Ri = first(_setdiff(clique, processed))
        push!(processed, Ri)
        Ri in non_ID_nodes && error("The propensity for $Ri is not identifiable. Estimation cannot be conducted.")
        append!(clique, propensitys[Ri].selection_set)
        bad = _intersect(propensitys[Ri].selection_set, non_ID_nodes)
        isempty(bad) || error("The propensity for $(join(bad, ", ")) needed for estimation is not identifiable.")
        clique = _unique(clique)
    end
    fit_rows = _all_one_rows(data, clique)
    sum(fit_rows) > 0 || error("No rows left for estimation. A larger data set is needed.")
    vocal && @info "$(sum(fit_rows)) rows are used in the M-estimation."
    ps = ones(Float64, nrow(data))
    for Rk in clique
        ps .*= [ismissing(x) ? 1.0 : Float64(x) for x in propensitys[Rk].pred]
    end
    weight = zeros(Float64, nrow(data))
    weight[fit_rows] .= 1 ./ ps[fit_rows]
    CliqueResult(clique, weight, fit_rows)
end

function M_estimation(data::DataFrame, est_FUN, cliques::Vector{CliqueResult}, start;
                      tol=1e-8, maxiter=100, fd_eps=1e-6)
    theta = Float64.(start)
    k = length(theta)
    n = nrow(data)

    function _raw_mat(th)
        raw = est_FUN(data, th)
        mat = raw isa AbstractVector ? reshape(Float64.(raw), :, 1) : Matrix{Float64}(raw)
        size(mat, 2) == length(cliques) || error("est_FUN must return one column per estimating equation/clique.")
        mat
    end

    function psi(th)
        mat = _raw_mat(th)
        vals = zeros(Float64, k)
        for j in eachindex(cliques)
            vals[j] = mean(mat[:, j] .* cliques[j].magic_weight)
        end
        vals
    end

    val = psi(theta)
    it = 0
    converged = norm(val) < tol
    while !converged && it < maxiter
        it += 1
        J = zeros(Float64, k, k)
        for j in 1:k
            step = zeros(Float64, k)
            step[j] = fd_eps * max(abs(theta[j]), 1.0)
            J[:, j] = (psi(theta .+ step) .- val) ./ step[j]
        end
        delta = J \ val
        theta .-= delta
        val = psi(theta)
        converged = norm(val) < tol || norm(delta) < tol
    end

    # Sandwich variance: Var(θ̂) = A⁻¹ B (A⁻¹)ᵀ / n
    # A = J (mean Jacobian of weighted estimating equations, already ≈ ∂ψ/∂θ)
    # B = (1/n) Σᵢ φᵢ φᵢᵀ  where φᵢ[j] = magic_weight_ij * est_FUN_ij
    J_final = zeros(Float64, k, k)
    val_final = psi(theta)
    for j in 1:k
        step = zeros(Float64, k)
        step[j] = fd_eps * max(abs(theta[j]), 1.0)
        J_final[:, j] = (psi(theta .+ step) .- val_final) ./ step[j]
    end
    phi = zeros(Float64, n, k)
    mat_final = _raw_mat(theta)
    for j in eachindex(cliques)
        phi[:, j] = mat_final[:, j] .* cliques[j].magic_weight
    end
    B_hat = phi' * phi / n
    vcov_theta, se_theta = try
        Jinv = inv(J_final)
        vc = Jinv * B_hat * Jinv' / n
        vc, sqrt.(max.(diag(vc), 0.0))
    catch
        fill(NaN, k, k), fill(NaN, k)
    end

    MEstimationResult(theta, vcov_theta, se_theta, converged, it, val_final)
end

function Mestimation(data::DataFrame, graph::MDAG, propensitys, ID::IDResult;
                     law=:target, est_FUN, variables, start, vocal=true, kwargs...)
    idlaw = _law(ID, law)
    nonid = [k for (k, v) in idlaw.id_status_collection if !v]
    length(variables) == length(start) || error("Length of variables and start must be the same.")
    cliques = CliqueResult[]
    for vars in variables
        if vars === nothing
            push!(cliques, CliqueResult(String[], ones(nrow(data)), trues(nrow(data))))
        else
            vset = vars isa AbstractString ? [String(vars)] : _unique(vars)
            estX = _intersect(vset, graph.missing_variables)
            estR = _findRX(graph.missing_indicators, graph.missing_variables, estX)
            push!(cliques, find_clique(graph, data, propensitys, estR, nonid; vocal=false))
        end
    end
    result = M_estimation(data, est_FUN, cliques, start; kwargs...)
    if vocal
        for i in eachindex(variables)
            se_str = isnan(result.se[i]) ? "SE=NA" : "SE=$(round(result.se[i], digits=4))"
            @info "$(sum(cliques[i].fit_rows)) rows used for parameter $i; estimate $(round(result.estimates[i], digits=3)), $se_str."
        end
    end
    (results=result, propensitys=propensitys, clique=cliques)
end

const f_ID_algorithm = ID_algorithm
const f_ID_flex_eval = ID_flex_eval
const f_propensity = propensity
const f_Mestimation = Mestimation
