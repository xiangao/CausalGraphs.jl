# CausalGraphs.jl — Project Notes

## Overview

Julia package for causal identification and semiparametric estimation on
ADMGs (acyclic directed mixed graphs — directed edges for causal relations,
bidirected edges for unobserved common causes), plus missing-data
identification/weighting on mDAGs. It is a from-scratch Julia rewrite of the
workflow in Anna Guo & Razieh Nabi's R packages `flexCausal` and
`flexMissing` (not a wrapper around them).

72/72 tests pass as of 2026-07-01 (`julia --project=. test/runtests.jl`,
~75s wall time). Testsets: graph (15), identification (12), finite-support ID
plug-in (12), optional NPCausal bridge (1), backdoor TMLE (5), NPS-TMLE (1),
nested ANIPW (6), missing data mDAG (15), missing weights (5).

## Architecture

**Routing pipeline** — the whole package is organized around one dispatch
point:

1. `make_graph()` builds an `ADMG` (`src/graph.jl`) — vertices/edges are
   `Symbol`s.
2. `identify(graph, treatment, outcome)` (`src/identification.jl`) checks, in
   order: a-fixable (`is_fix`, backdoor) → p-fixable (`is_p_fix`,
   front-door/NPS) → nested-fixable (`nested_fixability`, One-Line-ID style)
   → general Pearl-Shpitser `ID_algorithm` → `:not_identified`. Returns a
   named tuple whose `.strategy` field drives everything downstream.
3. `estimate_causal(...)` (`src/CausalGraphs.jl`) calls `identify` internally
   and routes to the matching estimator:
   - `:a_fixable` → `backdoor_tmle_a` (`src/backdoor.jl`) → TMLE, Onestep,
     Gcomp, IPW
   - `:p_fixable` → `nps_tmle_a` (`src/primalfix.jl`) → TMLE, Onestep
   - `:nested_fixable` → `nested_anipw_a` (`src/nested.jl`) → ANIPW, NIPW
   - `:id_algorithm` → `id_plugin_a` (`src/id_plugin.jl`) → finite-support
     plug-in with EIF-based confidence intervals (discrete variables only —
     it enumerates the observed support, it is **not** a general
     continuous-variable estimator)
4. Missing data is a separate, parallel pipeline: `make_mdag()` →
   `ID_algorithm(mdag)` (same exported name, different method, in
   `src/missing_core.jl`) → `propensity()` → `compute_missing_weights()`
   (`src/missing.jl`) → weights feed into `estimate_causal(...,
   sample_weights=wts)`.

**Key invariant — two separate graph types with different vertex types.**
`ADMG` (causal graphs, `src/graph.jl`) uses `Symbol` vertices/edges
(`Tuple{Symbol,Symbol}`). `MDAG` (missing-data graphs, `src/missing_core.jl`)
uses `String` vertices/edges (`Tuple{String,String}`). This is not a
Symbol-vs-String convenience wrapper — the two structs, and most of the
functions operating on them (`top_order`, `parents`, `children`,
`ID_algorithm`, ...), are genuinely separate implementations that happen to
share names via multiple dispatch. Don't assume a helper written for `ADMG`
works on `MDAG` or vice versa; check which struct a function's method
signature takes before reusing code across the two.

**`ID_algorithm` is overloaded across both worlds**: `ID_algorithm(graph,
treatment, outcome)` (ADMG, Pearl-Shpitser, symbolic) vs `ID_algorithm(mdag;
fulllaw=...)` (MDAG, missing-data target/full-law identification). Same
exported name, unrelated algorithms — check the argument types when reading
call sites.

**Treatments are binary, coded 0/1**, checked at estimator entry via
`check_binary_treatment!`. This is a hard current limitation, not just a
convention — nothing in the estimator code paths (backdoor, NPS, nested
ANIPW, ID plug-in) supports multi-valued or continuous treatment.

**`sample_weights` support is uneven.** Default/GLM nuisance fits accept
`sample_weights` throughout. Cross-fitted MLJ/SuperLearner nuisance fits
(`superlearner(binary=...)`, `src/ensemble.jl`) do not currently thread
weights through cross-fitting. If a session adds SuperLearner + weighted
estimation together, that gap needs to be closed or explicitly worked around.

## Relationship to CausalEstimate.jl and NPCausal.jl

- **CausalEstimate.jl** depends on `CausalGraphs` as a genuine weak
  dependency: its `Project.toml` has `[weakdeps] CausalGraphs = "68db26a2-..."`
  and `[extensions] CausalEstimateCausalGraphsExt = "CausalGraphs"`. This is
  wired correctly and is the "normal" direction of the relationship
  (CausalEstimate optionally builds on CausalGraphs' ADMGs/identification).

- **NPCausal.jl** also declares `CausalGraphs` as a weak dependency the same
  way (`[weakdeps]`/`[extensions]` in its own `Project.toml`), so that
  `NPCausal.admg_estimate_causal` can use `CausalGraphs`' identification.

- **The reverse bridge is broken as currently committed.** This repo ships
  `ext/CausalGraphsNPCausalExt.jl`, which defines
  `estimate_causal_npcausal_impl` and is meant to let
  `estimate_causal_npcausal(...)` (declared in `src/CausalGraphs.jl` via
  `Base.get_extension(@__MODULE__, :CausalGraphsNPCausalExt)`) delegate to
  `NPCausal.admg_estimate_causal`. **But this package's own `Project.toml`
  has no `[weakdeps]`/`[extensions]` section at all** — nothing declares
  `NPCausal` as a weak dep or wires the extension module. Verified directly:
  `Base.get_extension(CausalGraphs, :CausalGraphsNPCausalExt)` returns
  `nothing` even in a session where `NPCausal` is loaded, and `Pkg.develop`
  on NPCausal.jl silently added it as a normal `[deps]` entry rather than
  triggering extension resolution. Net effect: `estimate_causal_npcausal(...)`
  will **always** raise the "requires NPCausal.jl" error, regardless of
  whether NPCausal is installed/loaded, because the extension can never
  activate. Fix is to add to this package's `Project.toml`:
  ```toml
  [weakdeps]
  NPCausal = "cd612eeb-661c-4b56-8bb8-fdd3887649b5"

  [extensions]
  CausalGraphsNPCausalExt = "NPCausal"
  ```
  (UUID confirmed from NPCausal.jl's own `Project.toml`.) The
  `optional NPCausal bridge` testset in `test/runtests.jl` only exercises the
  error path (NPCausal not loaded) and would not have caught this — it never
  loads NPCausal and checks the success path.

- **UUID history**: this package's own UUID (`68db26a2-7c2b-457b-8f5d-
  f446513faf50`) was originally a placeholder and was fixed to a real UUID4 in
  commits `0bb1157` (update UUID) / `5ed7215` (fix self-UUID in Manifest),
  ~2026-06-14. Verified current: `CausalEstimate.jl/Project.toml`,
  `NPCausal.jl/Project.toml`, and this package's own `Project.toml` all now
  agree on `68db26a2-...`. Known failure mode already hit once: a stale
  `git-tree-sha1`/UUID pin in a downstream `docs/Manifest.toml` (fixed in
  CausalEstimate.jl's docs Manifest). If any *other* package in
  `~/projects/software/` fails to resolve/precompile `CausalGraphs` with a
  UUID mismatch error, check its `Manifest.toml` / `docs/Manifest.toml` for a
  leftover pin predating the UUID fix.

## Running tests

```bash
cd ~/projects/software/CausalGraphs.jl
julia --project=. test/runtests.jl
```

Single `test/runtests.jl`, ~75s, no `Pkg.test()` sandboxing issues observed
(unlike some sibling packages in this portfolio). Do not run
`Pkg.develop(path=...)` for another local package (e.g. NPCausal) inside this
repo's own `--project=.` environment to explore the extension bridge — it
mutates the tracked `Project.toml` by adding a normal `[deps]` entry instead
of resolving as a weak dep (this is itself the symptom of the missing
`[weakdeps]` wiring above). If you need to test the NPCausal bridge live, do
it in a scratch environment, not this repo's project.

## Docs

Documenter.jl site, deployed via `.github/workflows/docs.yml` on push to
`main` (also builds on PR without deploying). `docs/make.jl` lists 8
vignette pages plus `reference.md`:

- Graphs_and_Identification, ADMG_ID_Algorithm, Berkeley_Admissions_ID,
  Estimation_Backdoor_Frontdoor, Nested_and_Missing_Data,
  Smoking_Cessation_NHEFS, Job_Training_NSW, Job_Search_Mediation_JOBS

(README.md's vignette table currently lists 8 as well but names the fourth
one differently — "Real Example: Berkeley Admissions" etc. — cosmetic only,
same underlying files.) `deploydocs` targets `github.com/xiangao/
CausalGraphs.jl.git` with `devbranch="main"`. Building docs installs
Graphviz (`to_dot`/`draw_graph` DOT export) before `Pkg.develop(path=pwd())`.

## CI

`.github/workflows/CI.yml` currently runs the test matrix on **both**
`version: ["1.10", "1"]` (not just "1" — verify this hasn't changed
before assuming otherwise). `Project.toml` has `julia = "1.10"` in
`[compat]`, consistent with the matrix. If the matrix is ever trimmed to just
`"1"`, the `1.10` compat floor would go untested — check for that drift in
future sessions since this package sits at the bottom of a dependency chain
(CausalEstimate.jl, NPCausal.jl) that other work in this portfolio depends
on.

## Gotchas encountered

- `ADMG` is immutable; graph-transforming functions (`subgraph`,
  `fixed_graph`, `_remove_incoming` in identification.jl) return new `ADMG`
  values rather than mutating in place.
- `is_fix`/`is_p_fix`/nested-fixability/`ID_algorithm` are checked in a fixed
  priority order inside `identify()` — a graph that is nested-fixable is
  never reached if it also happens to be a-fixable or p-fixable, so
  `identify()`'s strategy reflects "cheapest applicable estimator," not
  "every valid identification strategy."
- Discrete/finite-support ID plug-in (`estimate_id`, `id_plugin_a`) is
  explicitly not a substitute for continuous-variable TMLE; it enumerates
  observed cells (see `IDPlugin_Y1.cell_EIF` test expecting exactly 8 cells
  for 3 binary variables in `test/runtests.jl`).
