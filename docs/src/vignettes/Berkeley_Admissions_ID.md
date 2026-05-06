# Berkeley Admissions with the ID Algorithm

```@meta
CurrentModule = CausalGraphs
```

This vignette summarizes a real-data finite-support ID example using the
Berkeley graduate admissions data. The full Quarto source is in
`vignettes/Berkeley_Admissions_ID.qmd`.

## Causal Question

What would the admission rate be if applicant gender were set to female rather
than male, under a front-door graph where department choice mediates the effect?

The data are the `UCBAdmissions` table distributed with R and mirrored by
[RDatasets](https://vincentarelbundock.github.io/Rdatasets/doc/datasets/UCBAdmissions.html).
They contain cell counts by admission outcome, gender, and department.

## Front-Door Graph

Let `female` be the treatment, `dept` the mediator, and `admitted` the outcome.
The applied front-door ADMG is:

```text
female -> dept -> admitted
female <-> admitted
```

```@example berkeley_id
using CausalGraphs

graph = make_graph(
    vertices = [:female, :dept, :admitted],
    di_edges = [(:female, :dept), (:dept, :admitted)],
    bi_edges = [(:female, :admitted)],
)

draw_graph(graph; direction="LR")
```

`identify()` routes this graph to the specialized front-door estimator class,
but the same graph is also identified by the general Pearl-Shpitser ID
algorithm:

```@example berkeley_id
route = identify(graph, :female, :admitted)
id = ID_algorithm(graph, :female, :admitted)

(strategy = route.strategy,
 identified_by_ID = id.identified,
 expression = string(id.expression))
```

## Data

The data are aggregated. `Freq` is used as a sample weight.

```@example berkeley_id
using DataFrames, DelimitedFiles, Downloads, Statistics

function read_rdatasets_csv(url)
    local_file = Downloads.download(url; timeout=120)
    x, header = readdlm(local_file, ',', Any, '\n'; header=true)
    DataFrame(x, Symbol.(vec(header)))
end

url = "https://vincentarelbundock.github.io/Rdatasets/csv/datasets/UCBAdmissions.csv"
raw = read_rdatasets_csv(url)

data = DataFrame(
    female = Float64.(raw.Gender .== "Female"),
    dept = String.(raw.Dept),
    admitted = Float64.(raw.Admit .== "Admitted"),
    weight = Float64.(raw.Freq),
)

(cells = nrow(data), applicants = Int(sum(data.weight)))
```

The raw weighted admission-rate contrast is:

```@example berkeley_id
weighted_mean(x, w) = sum(x .* w) / sum(w)

female_rate = weighted_mean(data.admitted[data.female .== 1],
                            data.weight[data.female .== 1])
male_rate = weighted_mean(data.admitted[data.female .== 0],
                          data.weight[data.female .== 0])

r = x -> round(x, sigdigits=4)
(female_rate = r(female_rate),
 male_rate = r(male_rate),
 raw_difference = r(female_rate - male_rate))
```

## ID Estimation

`estimate_id()` evaluates the symbolic ID functional by enumerating the observed
support and plugging in empirical conditional probabilities. It also reports a
finite-support EIF confidence interval. This is appropriate here because
`female`, `dept`, and `admitted` are finite-support variables.

```@example berkeley_id
id_res = estimate_id(
    a = [1.0, 0.0],
    data = select(data, Not(:weight)),
    graph = graph,
    treatment = :female,
    outcome = :admitted,
    sample_weights = data.weight,
)

(EYa1 = r(id_res[:IDPlugin_Y1].estimated_psi),
 EYa0 = r(id_res[:IDPlugin_Y0].estimated_psi),
 ACE = r(id_res[:IDPlugin].ACE),
 lower_ci = r(id_res[:IDPlugin].lower_ci),
 upper_ci = r(id_res[:IDPlugin].upper_ci),
 SE = r(id_res[:IDPlugin].standard_error),
 total_probability_a1 = r(id_res[:IDPlugin].total_probability_a1),
 total_probability_a0 = r(id_res[:IDPlugin].total_probability_a0))
```

The ID estimate answers the interventional question encoded by the front-door
graph. It is not the same object as the raw admission-rate contrast.

## Sensitivity Graph

If we add a direct `female -> admitted` edge while keeping the hidden
`female <-> admitted` association, the effect is no longer identified:

```@example berkeley_id
sensitivity_graph = make_graph(
    vertices = [:female, :dept, :admitted],
    di_edges = [(:female, :dept), (:dept, :admitted), (:female, :admitted)],
    bi_edges = [(:female, :admitted)],
)

identify(sensitivity_graph, :female, :admitted)
```

```@example berkeley_id
draw_graph(sensitivity_graph; direction="LR")
```

The point is the workflow: draw the graph, run `ID_algorithm()`, estimate only
when the graph identifies the effect, and then check how the conclusion changes
under a plausible alternative graph.
