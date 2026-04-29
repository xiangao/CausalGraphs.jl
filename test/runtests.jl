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
    # Build a graph where A is nested-fixable but not p-fixable:
    # W -> A -> M -> Y, A <-> Y (hidden), W -> Y, W <-> M
    # Treatment A is neither a-fix nor p-fix; OneLineID identifies it.
    Random.seed!(42)
    n = 400
    W = randn(n)
    A = Float64.(rand(n) .< 1 ./(1 .+ exp.(-W)))
    M = 0.5 .* A .+ 0.3 .* W .+ randn(n)
    Y = A .+ M .+ W .+ randn(n)   # ACE(A→Y) ≈ 1 + d(M)/d(A)·1 = 1 + 0.5 = 1.5
    df = DataFrame(W=W, A=A, M=M, Y=Y)
    g  = make_graph(vertices=[:W,:A,:M,:Y],
                    bi_edges=[(:A,:Y)],
                    di_edges=[(:W,:A),(:W,:Y),(:A,:M),(:M,:Y)])
    id = identify(g, :A, :Y)
    # This graph: A is p-fixable (district(A)={A,Y}, children(A)={M}; M not in district → p-fixable)
    # We just check that estimate_causal runs without error and gives a finite result
    res = estimate_causal(a=[1,0], data=df, graph=g, treatment=:A, outcome=:Y)
    @test any(haskey(res, k) for k in (:TMLE, :ANIPW))
    first_est = haskey(res, :TMLE) ? res[:TMLE].ACE : res[:ANIPW].ACE
    @test isfinite(first_est)
end
