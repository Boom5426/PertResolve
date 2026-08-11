# Deprecated results (do not use)

## `rep_grid_summary_5heads.csv` (moved here 2026-07-30)

The representation grid as originally committed. It is labelled a six-head result but was
computed on **five** heads. `rep_grid.py` built the MLP as
`MLPRegressor((128, 64), max_iter=500, random_state=0)`, and in scikit-learn >= 1.7 the first
positional parameter of `MLPRegressor` is `loss`, not `hidden_layer_sizes`. Every MLP fit
therefore raised `InvalidParameterError` and was discarded by a bare `except Exception:
continue`, so the neural-network family silently dropped out of all seven representations.
The same defect was present in `am_analysis.py`.

Both scripts now pass the parameter by keyword and count head-fit failures explicitly; the new
`reviewer_controls/rep_grid_summary.csv` carries `n_heads_scored`, `heads` and
`n_head_failures` columns so the failure cannot recur silently. Re-run on 2026-07-30 with all
six heads and zero failures:

| representation | PDS_cos 5 heads | PDS_cos 6 heads | Pearson-delta 5 heads | Pearson-delta 6 heads |
|---|---|---|---|---|
| theta | 0.475 | 0.473 | 0.632 | 0.632 |
| esm1v | 0.464 | 0.462 | 0.640 | 0.643 |
| esm2_global | 0.465 | 0.463 | 0.645 | 0.648 |
| esm2_globaldelta | 0.465 | 0.463 | 0.646 | 0.648 |
| esm2_sitedelta | 0.468 | 0.467 | 0.631 | 0.632 |
| esm2_window16 | 0.466 | 0.464 | 0.646 | 0.648 |
| esm2_sitedelta+theta | 0.465 | 0.464 | 0.627 | 0.628 |

The chance-level verdict is unchanged; every value moved by at most 0.003. The only reported
range that shifts is the PDS upper bound, from 0.475 (rounds to 0.48) to 0.473 (rounds to 0.47).

The main evaluation grid is **not** affected: `scripts/run_all_splits.py` and
`scripts/figures/all_splits_v4.py` construct the MLP with the keyword form, and the committed
`MLP-theta`, `MLP-esm` and `MLP-esm+theta` rows in the canonical tables carry real values.

## `second_probe_predictor_results.csv` (moved here 2026-07-10)

Superseded row-level LODO "rankability predictor" output. Its headline
`auroc` column (mean 0.974 core / 0.965 all) was shown by a multi-agent audit to be a
**config-identity + pseudo-replication artifact**, not a pilot-estimable per-perturbation
signal:

- The predictor was trained ROW-level over the ~18 correlated `(space x metric x n_work)`
  rows per perturbation (`n_test` equals the raw row count, e.g. Replogle 32517 rows over
  only 1832 perturbations). Its discrimination came from which distance metric / representation
  space a row used, not from the perturbation. Config-identity one-hots alone (zero biology)
  already reproduce TP53 0.955 / KRAS 0.936.
- The rankability label is exactly `rankable == (S > W)`, so `S, W, ratio, D_self, D_null`
  are the label; any predictor using them leaks the target.
- Two of the four features named in an earlier Methods draft ("PCA explained-variance ratio",
  "gene-space dimensionality") do not exist in the source table.
- At a genuine per-perturbation level TP53 and KRAS are uniformly un-rankable (single-class),
  so their committed AUROCs (0.943 / 0.964) are mathematically undefined and cannot exist honestly.

**Use instead:** `scripts/figures/rankability_predictor.py` ->
`results/rankability_predictor_honest.csv` (per-perturbation LODO, effect size as the sole
pilot-estimable feature, leaky/config features excluded, degenerate held-out datasets reported
as not evaluable). Honest result: mean AUROC 0.96 over the five non-degenerate datasets
(Replogle 0.97, Norman 0.93, Adamson 0.93, GATA1 0.95, JAK1 1.00); TP53 and KRAS not evaluable.
