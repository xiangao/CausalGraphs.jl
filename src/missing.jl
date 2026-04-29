# ── Missing data integration (v0.1 stub) ─────────────────────────────────────
#
# Full mDAG identification and propensity estimation is implemented in the
# standalone FlexMissing.jl package. This module provides a thin bridge:
# compute_missing_weights() runs FlexMissing's pipeline and returns per-row
# IPW weights that can be passed to the main estimation functions.
#
# v0.2 will embed FlexMissing's full graph engine here for joint estimation.

"""
    compute_missing_weights(mdag, data)

Given a FlexMissing.jl `MDAG` and a `DataFrame` with missing-indicator columns
(0/1), run the identification algorithm and estimate propensity scores, then
return a vector of inverse-probability weights (one per row).

Requires `FlexMissing` to be loaded in the calling session.

# Example
```julia
using FlexMissing, CausalGraphs, DataFrames
mdag = FlexMissing.make_graph(
    obs_variables=["A","Y"],
    missing_variables=["X"], missing_indicators=["Rx"],
    di_edges=[("A","Y"),("X","A"),("X","Y"),("X","Rx")]
)
ID  = FlexMissing.ID_algorithm(mdag)
wts = compute_missing_weights(mdag, data; ID=ID)
result = estimate_causal(a=[1,0], data=data, ..., sample_weights=wts)
```
"""
function compute_missing_weights(mdag, data::DataFrame; ID=nothing, kwargs...)
    error("""
Missing data integration requires FlexMissing.jl to be loaded and an `ID` result passed in.
Example:
    using FlexMissing
    ID  = FlexMissing.ID_algorithm(mdag)
    wts = compute_missing_weights(mdag, data; ID=ID)

Full embedded mDAG support is planned for CausalGraphs.jl v0.2.
""")
end
