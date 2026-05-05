# Lightweight graph visualization/export helpers.

struct MermaidGraph
    code::String
end

struct SVGGraph
    svg::String
end

"""
    to_dot(graph; direction="LR")

Return a Graphviz DOT representation of an ADMG.

Directed edges are blue. Bidirected edges, which represent hidden common
causes, are red and drawn with arrowheads at both ends. Fixed vertices are
drawn with square node shapes.
"""
function to_dot(g::ADMG; direction="LR")
    rankdir = _layout_direction(direction)
    lines = String["digraph {"]
    rankdir === nothing || push!(lines, "  graph [rankdir=\"$rankdir\"]")

    for v in g.vertices
        shape = get(g.fixed, v, false) ? "square" : "plaintext"
        push!(lines, "  $(_dot_id(v)) [shape=$shape, height=\".5\", width=\".5\"]")
    end

    for (from, to) in g.directed_edges
        push!(lines, "  $(_dot_id(from)) -> $(_dot_id(to)) [color=blue]")
    end
    for (from, to) in g.bidirected_edges
        push!(lines, "  $(_dot_id(from)) -> $(_dot_id(to)) [dir=both, color=red]")
    end

    push!(lines, "}")
    join(lines, "\n")
end

function to_dot(g::MDAG; direction="LR")
    rankdir = _layout_direction(direction)
    lines = String["digraph {"]
    rankdir === nothing || push!(lines, "  graph [rankdir=\"$rankdir\"]")

    for v in g.vertices
        shape = get(g.fixed, v, false) ? "square" : "plaintext"
        push!(lines, "  $(_dot_id(v)) [shape=$shape, height=\".5\", width=\".5\"]")
    end

    for (from, to) in g.di_edges
        push!(lines, "  $(_dot_id(from)) -> $(_dot_id(to)) [color=blue]")
    end

    push!(lines, "}")
    join(lines, "\n")
end

"""
    to_mermaid(graph; direction="LR")

Return a Mermaid flowchart representation of an ADMG.

This is convenient in Quarto, Markdown, and notebook contexts that support
Mermaid diagrams.
"""
function to_mermaid(g::ADMG; direction="LR")
    dir = _layout_direction(direction)
    dir === nothing && (dir = "LR")

    ids = Dict(v => _mermaid_id(v, i) for (i, v) in enumerate(g.vertices))
    lines = String["flowchart $dir"]

    for v in g.vertices
        push!(lines, "  $(ids[v])[\"$(_mermaid_label(v))\"]")
    end

    edge_styles = String[]
    edge_index = 0
    for (from, to) in g.directed_edges
        push!(lines, "  $(ids[from]) --> $(ids[to])")
        push!(edge_styles, "  linkStyle $edge_index stroke:#2563eb,stroke-width:1.8px")
        edge_index += 1
    end
    for (from, to) in g.bidirected_edges
        push!(lines, "  $(ids[from]) <--> $(ids[to])")
        push!(edge_styles, "  linkStyle $edge_index stroke:#dc2626,stroke-width:1.8px")
        edge_index += 1
    end

    push!(lines, "  classDef fixed stroke:#111827,stroke-width:2px")
    fixed_ids = [ids[v] for v in g.vertices if get(g.fixed, v, false)]
    isempty(fixed_ids) || push!(lines, "  class $(join(fixed_ids, ",")) fixed")

    append!(lines, edge_styles)
    join(lines, "\n")
end

function to_mermaid(g::MDAG; direction="LR")
    dir = _layout_direction(direction)
    dir === nothing && (dir = "LR")

    ids = Dict(v => _mermaid_id(v, i) for (i, v) in enumerate(g.vertices))
    lines = String["flowchart $dir"]

    for v in g.vertices
        push!(lines, "  $(ids[v])[\"$(_mermaid_label(v))\"]")
    end

    edge_styles = String[]
    for (edge_index, (from, to)) in enumerate(g.di_edges)
        push!(lines, "  $(ids[from]) --> $(ids[to])")
        push!(edge_styles, "  linkStyle $(edge_index - 1) stroke:#2563eb,stroke-width:1.8px")
    end

    push!(lines, "  classDef fixed stroke:#111827,stroke-width:2px")
    fixed_ids = [ids[v] for v in g.vertices if get(g.fixed, v, false)]
    isempty(fixed_ids) || push!(lines, "  class $(join(fixed_ids, ",")) fixed")

    append!(lines, edge_styles)
    join(lines, "\n")
end

"""
    to_svg(graph; direction="LR")

Render the graph to an SVG string using Graphviz (`dot`).
Returns `nothing` if `dot` is not available on the system PATH.
"""
function to_svg(g::Union{ADMG,MDAG}; direction="LR")
    dot_str = to_dot(g; direction=direction)
    try
        raw = read(pipeline(IOBuffer(dot_str), `dot -Tsvg`), String)
        # Strip XML declaration and DOCTYPE — keep only the <svg>...</svg> fragment
        m = match(r"(<svg\b.*</svg>)"s, raw)
        return m === nothing ? raw : m[1]
    catch
        return nothing
    end
end

"""
    draw_graph(graph; direction="LR")

Return a displayable graph object for an ADMG or MDAG.

Uses Graphviz (`dot`) to render an SVG when available — the SVG is embedded
directly in HTML with no JavaScript dependency. Falls back to a Mermaid
source object for environments without `dot` installed.
"""
function draw_graph(g::Union{ADMG,MDAG}; direction="LR")
    svg = to_svg(g; direction=direction)
    svg !== nothing && return SVGGraph(svg)
    return MermaidGraph(to_mermaid(g; direction=direction))
end

Base.show(io::IO, g::SVGGraph) = print(io, "[SVG graph — display in an HTML environment]")
Base.show(io::IO, ::MIME"text/plain", g::SVGGraph) = print(io, "[SVG graph — display in an HTML environment]")
Base.show(io::IO, ::MIME"text/html", g::SVGGraph) = print(io, g.svg)

Base.show(io::IO, g::MermaidGraph) = print(io, g.code)
Base.show(io::IO, ::MIME"text/plain", g::MermaidGraph) = print(io, g.code)

function Base.show(io::IO, ::MIME"text/html", g::MermaidGraph)
    print(io, "<pre class=\"mermaid\">\n")
    print(io, _escape_html(g.code))
    print(io, "\n</pre>")
end

function _layout_direction(direction)
    direction === nothing && return nothing
    dir = uppercase(String(direction))
    dir in ("LR", "RL", "TB", "BT", "TD") ||
        error("direction must be one of \"LR\", \"RL\", \"TB\", \"BT\", \"TD\", or nothing.")
    dir == "TD" ? "TB" : dir
end

_dot_id(v) = "\"" * replace(String(v), "\\" => "\\\\", "\"" => "\\\"") * "\""

function _mermaid_id(v, i::Int)
    raw = String(v)
    id = replace(raw, r"[^A-Za-z0-9_]" => "_")
    isempty(id) && (id = "node")
    occursin(r"^[A-Za-z_]", id) || (id = "v_" * id)
    "n$(i)_$id"
end

function _mermaid_label(v)
    replace(String(v), "\\" => "\\\\", "\"" => "#quot;", "[" => "(", "]" => ")")
end

function _escape_html(s::String)
    replace(s, "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")
end
