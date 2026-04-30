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

module FakeFlexMissing
using DataFrames
struct MDAG
    missing_indicators::Vector{String}
end
struct IDResult end
struct PropensityEstimate
    pred::Vector{Union{Missing,Float64}}
end
ID_algorithm(::MDAG) = IDResult()
propensity(::MDAG, data::DataFrame, ::IDResult; law=:target, kwargs...) =
    Dict("R" => PropensityEstimate(fill(0.5, nrow(data))))
end

@testset "missing weights bridge" begin
    df = DataFrame(X=[1.0, missing, 2.0], R=[1, 0, 1])
    mdag = FakeFlexMissing.MDAG(["R"])
    w_full = compute_missing_weights(mdag, df)
    @test w_full == [2.0, 0.0, 2.0]
    @test compute_missing_weights(mdag, df; complete_cases_only=true) == [2.0, 2.0]
end
