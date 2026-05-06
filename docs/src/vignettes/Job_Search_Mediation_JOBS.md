# Job Search Mediation JOBS II

```@meta
CurrentModule = CausalGraphs
```

This vignette uses the JOBS II job-search intervention data from the
`mediation` R package, mirrored by
[RDatasets](https://vincentarelbundock.github.io/Rdatasets/doc/mediation/jobs.html).
JOBS II was a randomized field experiment for unemployed workers. The
intervention offered job-skills workshops that taught job-search techniques and
coping strategies for setbacks in the job-search process. The follow-up data
record post-treatment depressive symptoms and employment-related measures. The
full Quarto source is in `vignettes/Job_Search_Mediation_JOBS.qmd`.

## Causal Question

How much of the effect of a job-skills workshop on depressive symptoms is
mediated by job-search self-efficacy?

The treatment is `treat`, the mediator is `job_seek`, the outcome is
`depress2`, and the baseline covariates are economic hardship, baseline
depression, sex, age, and nonwhite race.

## Data and Research Idea

The research idea is economic as well as behavioral. A job-search program might
affect mental health because it changes employment prospects, but it might also
work by changing how participants approach the job search itself. The mediator
`job_seek` captures job-search self-efficacy: whether participants feel more
capable of conducting a successful search.

The RDatasets version has 899 complete observations and is intended as an
illustrative analysis file. Here the goal is to show a complete causal
workflow: data, graph, identification, total-effect estimation, and a clearly
labeled mediation decomposition.

## Mediation DAG

Let `A` be workshop assignment, `M` be job-search self-efficacy, `Y` be
post-treatment depressive symptoms, and `X` be baseline variables.

```text
X -> A
X -> M
X -> Y
A -> M
M -> Y
A -> Y
```

```@example jobs_mediation
using CausalGraphs

conceptual_graph = make_graph(
    vertices = [:X, :A, :M, :Y],
    di_edges = [
        (:X, :A), (:X, :M), (:X, :Y),
        (:A, :M), (:M, :Y), (:A, :Y),
    ],
)

draw_graph(conceptual_graph; direction="LR")
```

For estimation, expand `X` into the observed baseline columns.

```@example jobs_mediation
baseline = [:econ_hard, :depress1, :sex, :age, :nonwhite]
vertices = [:treat, baseline..., :job_seek, :depress2]

edges = vcat(
    [(x, :treat) for x in baseline],
    [(x, :job_seek) for x in baseline],
    [(x, :depress2) for x in baseline],
    [(:treat, :job_seek), (:job_seek, :depress2), (:treat, :depress2)],
)

mediation_graph = make_graph(vertices=vertices, di_edges=edges)
identify(mediation_graph, :treat, :depress2).strategy
```

```@example jobs_mediation
draw_graph(mediation_graph; direction="LR")
```

```@example jobs_mediation
markov_pillow(mediation_graph, :treat; treatment=:treat)
```

For the total effect of `treat` on `depress2`, the mediator is part of the
causal pathway and is not adjusted away. The graph says that adjustment for the
baseline variables is sufficient.

## Total Effect Estimation

```@example jobs_mediation
using DataFrames, DelimitedFiles, Downloads, Statistics

function read_jobs()
    url = "https://vincentarelbundock.github.io/Rdatasets/csv/mediation/jobs.csv"
    local_file = Downloads.download(url; timeout=120)
    x, header = readdlm(local_file, ',', Any, '\n'; header=true)
    raw = DataFrame(x, Symbol.(vec(header)))

    DataFrame(
        treat = Float64.(raw.treat),
        econ_hard = Float64.(raw.econ_hard),
        depress1 = Float64.(raw.depress1),
        sex = Float64.(raw.sex),
        age = Float64.(raw.age),
        nonwhite = Float64.(raw.nonwhite .== "non.white1"),
        job_seek = Float64.(raw.job_seek),
        depress2 = Float64.(raw.depress2),
    )
end

data = read_jobs()

(n = nrow(data),
 treated = Int(sum(data.treat .== 1)),
 controls = Int(sum(data.treat .== 0)))
```

```@example jobs_mediation
total_res = estimate_causal(
    a = [1, 0],
    data = data,
    graph = mediation_graph,
    treatment = :treat,
    outcome = :depress2,
)

r = x -> round(x, sigdigits=4)

(TMLE_total_effect = r(total_res[:TMLE].ACE),
 lower_ci = r(total_res[:TMLE].lower_ci),
 upper_ci = r(total_res[:TMLE].upper_ci))
```

Lower values of `depress2` correspond to fewer depressive symptoms. Negative
effects therefore mean reductions in depressive symptoms. In this simple
analysis, the total-effect confidence interval crosses zero.

## Mediation Decomposition

For mediation, define `M0` as job-search self-efficacy under no workshop, `M1`
as job-search self-efficacy under the workshop, and `Y(a, m)` as depressive
symptoms under assignment `a` and mediator value `m`.

```text
total effect            = E[Y(1, M1) - Y(0, M0)]
natural direct effect   = E[Y(1, M0) - Y(0, M0)]
natural indirect effect = E[Y(1, M1) - Y(1, M0)]
```

The graph-based total effect above only requires identifying the intervention
`do(A=a)`. Natural direct and indirect effects require additional assumptions:
no unmeasured treatment-outcome, treatment-mediator, or mediator-outcome
confounding, plus the usual cross-world condition needed for natural effects.

The following code computes a parametric mediation g-formula:

```@example jobs_mediation
using LinearAlgebra

function mediator_matrix(d, avec)
    n = nrow(d)
    hcat(ones(n), avec, d.econ_hard, d.depress1, d.sex, d.age, d.nonwhite)
end

function outcome_matrix(d, avec, mvec)
    n = nrow(d)
    hcat(
        ones(n), avec, mvec, avec .* mvec,
        d.econ_hard, d.depress1, d.sex, d.age, d.nonwhite,
    )
end

n = nrow(data)
zeros_a = zeros(n)
ones_a = ones(n)

mediator_coef = mediator_matrix(data, data.treat) \ data.job_seek
m0 = mediator_matrix(data, zeros_a) * mediator_coef
m1 = mediator_matrix(data, ones_a) * mediator_coef

outcome_coef = outcome_matrix(data, data.treat, data.job_seek) \ data.depress2

ey_1m1 = mean(outcome_matrix(data, ones_a, m1) * outcome_coef)
ey_1m0 = mean(outcome_matrix(data, ones_a, m0) * outcome_coef)
ey_0m0 = mean(outcome_matrix(data, zeros_a, m0) * outcome_coef)

total_effect = ey_1m1 - ey_0m0
natural_direct = ey_1m0 - ey_0m0
natural_indirect = ey_1m1 - ey_1m0
prop_mediated = natural_indirect / total_effect

(ey_1m1 = r(ey_1m1),
 ey_1m0 = r(ey_1m0),
 ey_0m0 = r(ey_0m0),
 total_effect = r(total_effect),
 natural_direct = r(natural_direct),
 natural_indirect = r(natural_indirect),
 prop_mediated = r(prop_mediated))
```

In this linear plug-in analysis, most of the estimated reduction in depressive
symptoms is direct. About one fifth of the estimated reduction is mediated
through job-search self-efficacy.

## Optional Crumble.jl Handoff

This vignette does not need `Crumble.jl` to make the main point. The purpose
here is to show how `CausalGraphs.jl` moves from a DAG to total-effect
identification and estimation, and then to show the extra assumptions behind a
simple natural-effect mediation decomposition.

`Crumble.jl` is the better tool when the mediation estimand itself is the main
target. It implements mediation estimators for natural, organic, randomized
interventional, and recanting-twins effects, with nuisance estimation machinery
that is beyond the scope of this graph vignette. For the simple JOBS graph
above, the analogous `Crumble.jl` call would target natural effects:

```julia
using Crumble

crumble_res = crumble(
    data,
    ["treat"];
    outcome = "depress2",
    mediators = ["job_seek"],
    covar = string.(baseline),
    effect = "N",
)
```

This is intentionally shown as optional code rather than a required dependency
of `CausalGraphs.jl`. The graph package should be usable without pulling in a
full mediation-estimation stack. In a substantive mediation analysis, however,
one natural workflow is to use `CausalGraphs.jl` to state and stress-test the
DAG assumptions, use `identify()` to separate total-effect identification from
mechanism assumptions, and then use `Crumble.jl` when the target is a formal
mediation estimand and the data support the required assumptions.

The distinction matters in the next section. If the mediator-outcome
relationship is confounded by an unmeasured variable, the total effect may
still be identifiable from the graph, but the natural direct and indirect
effects are no longer justified by the simple natural-effect formula. A
different mediation estimand or additional measured mediator-outcome
confounders would be needed before using a richer mediation estimator.

## Mediator-Outcome Confounding

The weakest link in many mediation analyses is the mediator-outcome
relationship. Job-search self-efficacy is post-treatment, not randomized.
Unmeasured optimism, family support, health, or local job-market conditions
could affect both self-efficacy and depressive symptoms.

A sensitivity graph adds a bidirected edge between the mediator and outcome.
This is different from the directed edge `job_seek -> depress2`, which is the
causal mediated pathway. The bidirected edge `job_seek <-> depress2` means
there is an unmeasured common cause of job-search self-efficacy and depressive
symptoms.

```@example jobs_mediation
sensitivity_graph = make_graph(
    vertices = vertices,
    di_edges = edges,
    bi_edges = [(:job_seek, :depress2)],
)

draw_graph(sensitivity_graph; direction="LR")
```

The total effect of `treat` on `depress2` is still a treatment-effect question
and can still be estimated under the baseline-adjustment assumptions above,
because the total effect does not condition on or intervene on the mediator.
`CausalGraphs.jl` therefore still identifies the total effect in this
sensitivity graph:

```@example jobs_mediation
identify(sensitivity_graph, :treat, :depress2).strategy
```

```@example jobs_mediation
markov_pillow(sensitivity_graph, :treat; treatment=:treat)
```

```@example jobs_mediation
sensitivity_total_res = estimate_causal(
    a = [1, 0],
    data = data,
    graph = sensitivity_graph,
    treatment = :treat,
    outcome = :depress2,
)

(TMLE_total_effect = r(sensitivity_total_res[:TMLE].ACE),
 lower_ci = r(sensitivity_total_res[:TMLE].lower_ci),
 upper_ci = r(sensitivity_total_res[:TMLE].upper_ci))
```

The natural direct and indirect effects, however, are no longer justified by
the simple mediation g-formula because the `M -> Y` relationship is confounded
after conditioning on `A` and `X`.

## Takeaway

This example separates three tasks that are often collapsed in applied
mediation work:

- the DAG encodes the economic and behavioral story;
- `identify()` and `estimate_causal()` handle the total effect of the workshop;
- the natural-effect decomposition requires extra assumptions and an explicit
  mediation formula;
- `Crumble.jl` is the natural companion when the mediation estimand, rather
  than the graph-to-total-effect workflow, is the main analysis target.
