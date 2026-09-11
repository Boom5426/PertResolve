# Resolution panel v2

The committed panel artifacts below are frozen historical outputs. Their
`detectable`, `identifiable`, `ceiling` and `verdict` columns are legacy display
fields, not the current public API. Read `ceiling` as the empirical split-half
reproducibility reference, independently of detection and identification; it
is not a hard ceiling or bound, and it does not alone certify model-ranking
resolution. New aggregation uses the independent axis names documented in
[`results/canonical/README.md`](../canonical/README.md).

Frozen protocol: `docs/PREREG_RESOLUTION_PANEL_v2.md`. Read it before reading any number
here. Runner: `scripts/analysis/resolution_panel_v2.py`. Compute is remote; these files are
the synced outputs.

One dataset is one (resource, stratum) pair. Each run emits `<name>.json` (scalars plus the
provenance of that run) and `<name>_units.csv` (per-perturbation detection, identification
and their conjunction, which is what makes the jointly evaluable fraction derivable).

## Gate 1, loader reproduction: PASSED 2026-09-01

Section 8b requires the new loader to reproduce the v1 panel before anything else runs. It
reproduces it **exactly**, not to within tolerance, on every scalar of both datasets:

| dataset | n | detectable | identifiable | split-half reference | rho2 | rho2_nn_median | p_correct_order |
|---|---|---|---|---|---|---|---|
| norman2019 | 195 | 0.6871794871794872 | 0.005128205128205128 | 0.8975713719270421 | 5.947260807581858 | 0.7108467641514068 | 0.997 |
| adamson2016 | 94 | 0.4787234042553192 | 0.031914893617021274 | 0.8509637382749943 | 4.657357310246475 | 0.5481471366042009 | 0.974 |

Every value is bit-identical to the corresponding row of `results/canonical/resolution_panel.csv`.
That also demonstrates the section 7a change is inert: storing the per-perturbation
identification vector moved no reported quantity.

New in v2 and not available in v1: `jointly_evaluable_fraction`, 0.005128205128205128 for
Norman and 0.031914893617021274 for Adamson. In these two rows, it equals the historical
identifiable fraction, so every identified perturbation in this particular comparison also
passes the detection criterion. The equality is an observed intersection, not a logical rule;
the detection and identification axes remain independent elsewhere in the panel.

### One defect the gate caught, which is why it exists

The first version of the runner filtered out the cells of sub-threshold perturbations before
reducing dimension. That fitted the reduction on 106,022 of Norman's 111,445 cells instead of
all of them, and emptied the `excluded` provenance. It moved the empirical reference by 5.6e-5 and
`rho2_nn_median` by 1.6e-3 while leaving both discrete fractions bit-exact, so it would have
been easy to wave through as numerical noise. With no cap active the runner now hands the
criterion every cell, exactly as `pertresolve.resolution.cli` does, and lets
`group_profiles` do the excluding and report it.

## Runner provenance, and a schema drift declared rather than repaired

Recorded 2026-09-03, when this directory was brought under version control. **The runner now
in the repository is a later generation than the results here, and no re-run was permitted, so
the difference is stated instead of removed.**

What was checked, and how:

- The current `scripts/analysis/resolution_panel_v2.py` writes `source` and
  `extract_provenance` into each summary. **Zero of the 29 JSON files carry either key**, and
  all 29 carry `h5ad`, which the current runner no longer writes under that name. The runner on
  disk therefore did not produce these files.
- At least two generations are visible among the results themselves. `gate1/adamson2016.json`
  and `gate1/norman2019.json` carry a scalar `perturbation_key` and none of
  `reduce_path`, `n_cells_unassigned_dropped`, `unit_separator`, `dense_working_set_gb`. The
  other 27 carry all four and a list-valued `perturbation_key`. With the runner that is a
  third.

What this does and does not affect:

- It does **not** put the numbers in doubt. Gate 1 reproduces the v1 panel bit-exactly on every
  scalar, the section 7a change to `pertresolve/resolution/report.py` is provably inert
  (`identifiable_fraction` is still `identified.mean()`), and every value in
  `results/canonical/resolution_panel_v2_table.csv` is asserted equal to its source JSON by
  `scripts/analysis/build_resolution_panel_table.py` on every run.
- It **does** mean a reader cannot regenerate these files from this repository and expect the
  same JSON schema. Until the panel is re-run under the committed runner, treat the JSON files
  as the artifact and the runner as the current implementation of the same protocol, not as its
  reproduction recipe.

A re-run under the committed runner is the repair. It is recorded as a future upgrade in
`docs/FIG456_NARRATIVE_REDESIGN_FEASIBILITY.md` and was not executed.

## The guard that was missing

The runner accepted `--out` without calling `pertresolve.paths.reject_repo_results`, and it
wrote into `results/`, which that function exists to protect. A tracked test,
`tests/test_repo_contracts.py::test_every_out_taking_script_guards_the_results_tree`, was
failing because of it. The guard was added on 2026-09-03 and the test now passes. The outputs
already in this directory predate the guard; they were not moved, because moving them would
break the paths `RESULTS.md` cites.

## Measured cost, on an idle node (24 cores, 124 GB)

| dataset | cells | load | reduce | score | total |
|---|---|---|---|---|---|
| norman2019 | 111,445 | 4.7 s | 54.6 s | 7.7 s | 67.0 s |
| adamson2016 | 65,337 | 6.8 s | 43.8 s | 2.1 s | 52.7 s |

Dimension reduction dominates, which is what makes the section 13c candidate-pool cap the
operative compute lever on the three largest resources rather than the per-perturbation cell
cap.
