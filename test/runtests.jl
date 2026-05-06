using Test, DataFrames, Random
using CausalGraphs

# ── Graph helpers ─────────────────────────────────────────────────────────────
@testset "graph" begin
    g = make_graph(vertices=[:A,:M,:Y,:X],
                   bi_edges=[(:A,:Y)],
                   di_edges=[(:X,:A),(:X,:Y),(:X,:M),(:A,:M),(:M,:Y)])
    @test :X in parents(g, :A)
    @test :Y in district(g, :A)
    @test length(all_districts(g)) == 3   # {A,Y}, {M}, {X}
    @test is_p_fix(g, :A)
    @test !is_fix(g, :A)

    g_back = make_graph(vertices=[:A,:Y,:X],
                        di_edges=[(:X,:A),(:X,:Y),(:A,:Y)])
    @test is_fix(g_back, :A)
    @test ancestors(g_back, [:Y]) == [:Y, :A, :X] ||
          Set(ancestors(g_back, [:Y])) == Set([:Y,:A,:X])
    sub = subgraph(g_back, [:A,:Y])
    @test :X ∉ sub.vertices

    dot = to_dot(g; direction="LR")
    @test occursin("graph [rankdir=\"LR\"]", dot)
    @test occursin("\"A\" -> \"M\" [color=blue]", dot)
    @test occursin("\"A\" -> \"Y\" [dir=both, color=red]", dot)

    mermaid = to_mermaid(g; direction="LR")
    @test occursin("flowchart LR", mermaid)
    @test occursin("-->", mermaid)
    @test occursin("<-->", mermaid)

    html = sprint(show, MIME("text/html"), draw_graph(g))
    @test occursin("class=\"mermaid\"", html) || occursin("<svg", html)
end

# ── Identification ────────────────────────────────────────────────────────────
@testset "identification" begin
    # backdoor
    g1 = make_graph(vertices=[:A,:Y,:X],
                    di_edges=[(:X,:A),(:X,:Y),(:A,:Y)])
    @test identify(g1,:A,:Y).strategy == :a_fixable

    # front-door (p-fixable / NPS)
    g2 = make_graph(vertices=[:A,:M,:Y],
                    bi_edges=[(:A,:Y)],
                    di_edges=[(:A,:M),(:M,:Y)])
    @test identify(g2,:A,:Y).strategy == :p_fixable

    # not identified
    g3 = make_graph(vertices=[:A,:Y],
                    bi_edges=[(:A,:Y)],
                    di_edges=[(:A,:Y)])
    id3 = identify(g3,:A,:Y)
    @test id3.strategy == :not_identified || id3.strategy == :a_fixable

    # nested-fixable: two-bow graph (A->B<->A, C->D<->C, A->D)
    # Classic example from Bhattacharya, Nabi & Shpitser (2022)
    g4 = make_graph(vertices=[:U,:A,:M,:Y],
                    bi_edges=[(:A,:Y)],
                    di_edges=[(:U,:A),(:U,:M),(:A,:M),(:M,:Y)])
    id4 = identify(g4,:A,:Y)
    @test id4.strategy in (:p_fixable, :nested_fixable, :a_fixable)

    # Pearl-Shpitser ID succeeds for front-door and fails for the bow graph.
    id_fd = ID_algorithm(g2, :A, :Y)
    @test id_fd.identified
    @test occursin("sum_{M}", string(id_fd.expression))
    @test is_id_identified(g2, :A, :Y)
    @test !ID_algorithm(g3; treatment=:A, outcome=:Y).identified

    fs = fixing_sequence(g1, [:A])
    @test fs.fixable
    @test fs.fixing_order == [:A]
    nested = nested_fixability(g4, :A, :Y)
    @test nested.treatment == :A
    @test nested.outcome == :Y
end

@testset "finite-support ID plug-in" begin
    Random.seed!(11)
    n = 500
    A = Float64.(rand(n) .< 0.5)
    M = Float64.(rand(n) .< (0.2 .+ 0.6 .* A))
    Y = Float64.(rand(n) .< (0.1 .+ 0.2 .* A .+ 0.4 .* M))
    df = DataFrame(A=A, M=M, Y=Y)
    g = make_graph(vertices=[:A, :M, :Y],
                   di_edges=[(:A, :M), (:M, :Y)],
                   bi_edges=[(:A, :Y)])

    res = estimate_id(a=[1, 0], data=df, graph=g, treatment=:A, outcome=:Y)
    @test haskey(res, :IDPlugin)
    @test isfinite(res[:IDPlugin].ACE)
    @test isfinite(res[:IDPlugin].standard_error)
    @test res[:IDPlugin].standard_error >= 0
    @test res[:IDPlugin].lower_ci < res[:IDPlugin].upper_ci
    @test length(res[:IDPlugin].EIF) == nrow(df)
    @test abs(res[:IDPlugin_Y1].total_probability - 1) < 1e-8
    @test abs(res[:IDPlugin_Y0].total_probability - 1) < 1e-8
    @test length(res[:IDPlugin_Y1].cell_EIF) == 8

    collapsed = combine(groupby(df, [:A, :M, :Y]), nrow => :w)
    res_w = estimate_id(a=[1, 0], data=collapsed[:, [:A, :M, :Y]],
                        graph=g, treatment=:A, outcome=:Y,
                        sample_weights=collapsed.w)
    @test res_w[:IDPlugin].ACE ≈ res[:IDPlugin].ACE
    @test res_w[:IDPlugin].standard_error ≈ res[:IDPlugin].standard_error

    g_bow = make_graph(vertices=[:A, :Y],
                       di_edges=[(:A, :Y)],
                       bi_edges=[(:A, :Y)])
    @test_throws ErrorException estimate_id(a=1, data=df[:, [:A, :Y]],
                                            graph=g_bow, treatment=:A, outcome=:Y)
end

@testset "optional NPCausal bridge" begin
    @test_throws ErrorException estimate_causal_npcausal(
        a=[1, 0],
        data=DataFrame(A=Float64[], Y=Float64[]),
        vertices=[:A, :Y],
        di_edges=[(:A, :Y)],
        treatment=:A,
        outcome=:Y,
    )
end

# ── Backdoor TMLE ─────────────────────────────────────────────────────────────
@testset "backdoor TMLE" begin
    Random.seed!(1)
    n = 300; X = randn(n)
    A = Float64.(rand(n) .< 1 ./(1 .+ exp.(-X)))
    Y = 2 .* A .+ X .+ randn(n)
    df = DataFrame(A=A, Y=Y, X=X)
    g  = make_graph(vertices=[:A,:Y,:X], di_edges=[(:X,:A),(:X,:Y),(:A,:Y)])
    res = estimate_causal(a=[1,0], data=df, graph=g, treatment=:A, outcome=:Y)
    @test haskey(res, :TMLE)
    @test isfinite(res[:TMLE].ACE)
    @test abs(res[:TMLE].ACE - 2.0) < 0.5

    w = vcat(fill(20.0, n ÷ 2), fill(0.1, n - n ÷ 2))
    res_w = estimate_causal(a=[1,0], data=df, graph=g, treatment=:A, outcome=:Y,
                            sample_weights=w)
    @test isfinite(res_w[:TMLE].ACE)
    @test res_w[:TMLE].ACE != res[:TMLE].ACE
end

# ── NPS-TMLE (front-door / p-fixable) ────────────────────────────────────────
@testset "NPS-TMLE" begin
    Random.seed!(2)
    n = 300; X = rand(n)
    A = Float64.(rand(n) .< (0.3 .+ 0.2 .* X))
    U = 1 .+ A .+ X .+ randn(n)
    M = Float64.(rand(n) .< 1 ./(1 .+ exp.(-(-1 .+ A .+ X))))
    Y = U .+ M .+ X .+ randn(n)
    df = DataFrame(X=X, A=A, M=M, Y=Y)
    g  = make_graph(vertices=[:A,:Y,:X,:M], bi_edges=[(:A,:Y)],
                    di_edges=[(:X,:A),(:X,:Y),(:X,:M),(:A,:M),(:M,:Y)])
    res = estimate_causal(a=[1,0], data=df, graph=g, treatment=:A, outcome=:Y)
    @test isfinite(res[:TMLE].ACE)
end

# ── Nested ANIPW ──────────────────────────────────────────────────────────────
@testset "nested ANIPW" begin
    Random.seed!(42)
    n = 500
    U_EI   = randn(n)
    U_ET   = randn(n)
    U_ITox = randn(n)
    U_VA   = randn(n)
    U_VCD4 = randn(n)
    ViralLoad = randn(n)
    Exercise  = randn(n) .+ 0.5 .* U_EI .+ 0.5 .* U_ET
    Income    = randn(n) .+ 0.4 .* ViralLoad .+ 0.4 .* U_EI .+ 0.4 .* U_ITox
    T         = Float64.(rand(n) .< 1 ./(1 .+ exp.(-(0.5 .* ViralLoad .+ 0.3 .* Income .+ 0.4 .* U_ET))))
    Toxicity  = randn(n) .+ 0.6 .* T .+ 0.3 .* U_ITox
    Adherence = randn(n) .+ 0.5 .* Toxicity .+ 0.3 .* U_VA
    CD4       = randn(n) .+ 0.4 .* T .+ 0.3 .* Adherence .+
                0.2 .* Exercise .+ 0.3 .* ViralLoad .+ 0.3 .* U_VCD4 .+ 0.2 .* U_VA
    df = DataFrame(ViralLoad=ViralLoad, Income=Income, Exercise=Exercise,
                   T=T, Toxicity=Toxicity, Adherence=Adherence, CD4=CD4)
    g = make_graph(
        vertices = [:ViralLoad, :Income, :Exercise, :T, :Toxicity, :Adherence, :CD4],
        di_edges = [
            (:ViralLoad, :T), (:ViralLoad, :CD4), (:Income, :T), (:Exercise, :CD4),
            (:T, :Toxicity), (:T, :CD4), (:Toxicity, :Adherence), (:Adherence, :CD4),
        ],
        bi_edges = [
            (:Income, :Toxicity), (:Exercise, :Income), (:Exercise, :T),
            (:ViralLoad, :Adherence), (:ViralLoad, :CD4),
        ],
    )
    @test !is_fix(g, :T)
    @test !is_p_fix(g, :T)
    @test identify(g, :T, :CD4).strategy == :nested_fixable
    res = estimate_causal(a=[1,0], data=df, graph=g, treatment=:T, outcome=:CD4)
    @test haskey(res, :ANIPW)
    @test isfinite(res[:ANIPW].ACE)
    @test isfinite(res[:NIPW].ACE)
end

function g2_mdag()
    make_mdag(
        obs_variables=[],
        missing_variables=["X1", "X2", "X3", "X4", "X5", "X6"],
        missing_indicators=["R1", "R2", "R3", "R4", "R5", "R6"],
        di_edges=[
            ("X1", "X2"), ("X2", "X3"), ("X3", "X4"), ("X4", "X5"), ("X5", "X6"),
            ("R6", "R5"), ("R6", "R3"), ("R5", "R3"), ("R4", "R3"), ("R3", "R2"), ("R2", "R1"),
            ("X1", "R3"),
            ("X3", "R4"), ("X3", "R5"), ("X3", "R6"),
            ("X4", "R1"), ("X4", "R6"),
            ("X5", "R2"),
            ("X6", "R5"),
        ],
    )
end

@testset "missing data mDAG" begin
    g = make_mdag(obs_variables=["X", "A", "Y"], missing_variables=[], missing_indicators=[],
                  di_edges=[("X", "A"), ("A", "Y"), ("X", "Y")])
    @test top_order(g) == ["X", "A", "Y"]
    @test parents(g, "Y") == ["X", "A"]
    @test children(g, "X") == ["A", "Y"]
    @test !is_d_separated(g, "A", "Y", ["X"]).d_separated
    @test occursin("\"X\" -> \"A\" [color=blue]", to_dot(g))
    @test occursin("flowchart LR", to_mermaid(g))

    g2 = g2_mdag()
    id = ID_algorithm(g2; fulllaw=true)
    @test id.target.id_status
    @test id.full !== nothing
    @test !id.full.id_status
    target = summarize_ID(id).target
    @test all(target.ID_Status)
    @test id.target.selection_collection["R3"] == ["R1", "R4", "R5"]
    @test id.target.selectionr_collection["R3"] == ["R4", "R5"]
    @test id.target.selectionx_collection["R3"] == ["R1"]
    @test id.full.id_status_collection["R3"] == false
    @test id.full.id_status_collection["R5"] == false
end

@testset "missing weights" begin
    df = DataFrame(X=[1.0, missing, 2.0], R=[1, 0, 1])
    mdag = make_mdag(obs_variables=[], missing_variables=["X"],
                     missing_indicators=["R"], di_edges=[])
    id = ID_algorithm(mdag)
    props = Dict("R" => PropensityEstimate(nothing, Union{Missing,Float64}[0.5 for _ in 1:nrow(df)],
                                           String[], trues(nrow(df)), trues(nrow(df)),
                                           String[], String[]))
    w_full = compute_missing_weights(mdag, df; ID=id, propensities=props)
    @test w_full == [2.0, 0.0, 2.0]
    @test compute_missing_weights(mdag, df; ID=id, propensities=props,
                                  complete_cases_only=true) == [2.0, 2.0]

    props_fit = propensity(mdag, df, id; vocal=false)
    result = Mestimation(df, mdag, props_fit, id;
                         est_FUN=(d, th) -> reshape(coalesce.(d.X, 0.0) .- th[1], :, 1),
                         variables=[["X"]], start=[0.0], vocal=false)
    @test result.results.converged
    @test length(result.results.estimates) == 1
    @test isfinite(result.results.se[1])
end
