# ── Missing data integration ─────────────────────────────────────────────────

"""
    compute_missing_weights(mdag, data; ID=nothing, law=:target,
                            indicators=nothing, complete_cases_only=false,
                            normalize=false, kwargs...)

Given a CausalGraphs.jl `MDAG` and a `DataFrame` with missing-indicator columns
(0/1), run the missing-data identification algorithm and estimate propensity
scores, then return inverse-probability weights for complete-case estimation.

By default the returned vector has one entry per row of `data`, with zero
weight on incomplete rows. If `complete_cases_only=true`, only weights for rows
where all requested indicators equal 1 are returned; this is useful when passing
`dropmissing(data)` to `estimate_causal`.

# Example
```julia
using CausalGraphs, DataFrames
mdag = make_mdag(
    obs_variables=["A","Y"],
    missing_variables=["X"], missing_indicators=["Rx"],
    di_edges=[("A","Y"),("X","A"),("X","Y"),("X","Rx")]
)
ID  = ID_algorithm(mdag)
wts = compute_missing_weights(mdag, data; ID=ID, complete_cases_only=true)
result = estimate_causal(a=[1,0], data=dropmissing(data), ..., sample_weights=wts)
```
"""
function compute_missing_weights(mdag::MDAG, data::DataFrame; ID=nothing, law=:target,
                                 indicators=nothing, propensities=nothing,
                                 complete_cases_only=false, normalize=false,
                                 eps_prob=1e-10, kwargs...)
    ID = ID === nothing ? ID_algorithm(mdag) : ID
    props = propensities === nothing ?
            propensity(mdag, data, ID; law=law, kwargs...) :
            propensities

    Rs = indicators === nothing ? String.(getproperty(mdag, :missing_indicators)) :
         String.(indicators isa Union{AbstractString,Symbol} ? [indicators] : indicators)
    n = nrow(data)
    isempty(Rs) && return ones(Float64, n)

    missing_props = setdiff(Rs, String.(collect(keys(props))))
    isempty(missing_props) ||
        error("No missing-data propensity estimate was returned for $(join(missing_props, ", ")).")

    complete = trues(n)
    for R in Rs
        R in names(data) || error("Missing indicator column `$R` is not present in `data`.")
        complete .&= [!ismissing(v) && Float64(v) == 1.0 for v in data[!, R]]
    end

    weights = zeros(Float64, n)
    weights[complete] .= 1.0
    for R in Rs
        pred = props[R].pred
        length(pred) == n || error("Propensity predictions for `$R` do not match `data` rows.")
        for i in findall(complete)
            p = pred[i]
            ismissing(p) && error("Propensity prediction for `$R` is missing on complete row $i.")
            weights[i] /= max(Float64(p), eps_prob)
        end
    end

    if normalize
        denom = mean(weights[complete])
        denom > 0 && (weights[complete] ./= denom)
    end
    complete_cases_only ? weights[complete] : weights
end
