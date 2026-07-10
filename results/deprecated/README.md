# Deprecated results (do not use)

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
