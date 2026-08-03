# Corrected measurement-resolution scaling, 2026-08-03

Backs `results/canonical/floor_law_v2.csv`, produced by
`scripts/analysis/resolution_scaling.py` with the estimator in
`alleleperturb/resolution/scaling.py`. Supersedes `floor_law.csv` for the analytical
layer; `floor_law.csv` and its generator `floor_law_fit.py` are unchanged and still
reproduce, so the previously reported values remain traceable.

## Why the previous table could not support the claim it was cited for

Methods described a dimensionless discrimination signal-to-noise `rho` and asserted that
discrimination depends on the data only through it. `floor_law.csv` reported
`delta2 = 0, rho = 0` in 13 of 16 gene-by-depth rows, so the axis was a single point for
three of the four datasets and the assertion had three informative points, all from JAK1.
Two defects in the estimator caused this, and both are now fixed.

**The noise term was inflated by the wild-type mean.** The old estimator took

    eta2_hat = 0.5 * mean_v ||(B_v - w_B) - (E_v - w_E)||^2

with two *independent* wild-type half-means, while `pair2` was computed among profiles
that all share `w_E`, where it cancels exactly. The wild-type sampling noise therefore
entered the noise term twice and the signal term not at all, biasing `delta2` down by
about twice that noise. Measured on the real splits, `eta2_wt` is as large as the variant
sampling noise itself:

| gene | m | eta2_wt | eta2 as-is | eta2 corrected | delta2 as-is | delta2 corrected |
|---|---|---|---|---|---|---|
| TP53 | 25 | 36.90 | 74.59 | 37.63 | -73.92 | 0.02 |
| GATA1 | 25 | 103.90 | 201.94 | 98.57 | -203.42 | 3.33 |
| JAK1 | 25 | 72.44 | 151.46 | 79.15 | -60.80 | 83.84 |

so the estimate did not land slightly below zero, it landed near -74 and -203.

**The estimate was clipped at zero.** `delta2 = max(pair2 - 2 * eta2, 0.0)` then turned
those into exact zeros. On data with no true signal, half of all unbiased estimates are
negative by construction, so clipping reports certainty the data does not carry. Nothing
is clipped now; the sign is kept and an interval is reported with it.

## What the estimator does now

Each variant's cells are split into four disjoint groups of `m` cells: Q and T score the
replicate ceiling, S1 and S2 estimate the signal. Nothing is shared between the two axes,
which removes the shared noise realisation the old script had between them, at the cost of
requiring 4m rather than 2m cells per variant. That cost is visible as `n_var`.

* `delta2 = pair2 - 2 * eta2 / H`, unbiased and unclipped, with a **two-level bootstrap**
  over variants and over cell splits. Resampling variants alone gave intervals that
  excluded zero for a quantity that cannot be negative (KRAS fixed cohort, m = 25:
  [-0.511, -0.014]); with splits resampled the same row reads [-0.701, 0.191].
* A **permutation P value** for the null that no between-variant signal exists. A REML
  variance-components fit was considered and rejected: it constrains the between variance
  to be non-negative, reintroducing exactly the boundary artefact the second defect is
  about, so it cannot serve as an independent check. This substitutes for the REML step in
  the approved plan.
* The wild-type reference is built from whatever wild-type cells exist rather than
  silently dropping to no subtraction. TP53 carries fewer than 200 wild-type cells, so a
  rule demanding 2m of them left the m = 100 and m = 150 profiles as raw means, which a
  cosine score does not treat the same. `wt_cells_per_ref` and `wt_ref_full_depth` record
  what was actually used.

The estimator is validated against data with a known signal in
`tests/test_resolution_scaling.py`: unbiasedness across four signal regimes including
exactly zero, the reproduction of the wild-type bias (independent references recover 0.99
against a true 2.00, matching the predicted deficit of 1.00), permutation P uniform under
the null with 100% power at signal, and bootstrap coverage of 94.5% one-level and 96.7%
two-level against a nominal 95%.

## Result

The axis is no longer degenerate. Over the 16 eligible-cohort points:

| gene | rho2 range | permutation P | ceiling |
|---|---|---|---|
| TP53 | -0.005 to -0.001, all intervals covering 0 | 0.58 to 0.84 | 0.473 to 0.499 |
| KRAS | -0.002 to 0.004, all intervals covering 0 | 0.14 to 0.76 | 0.495 to 0.510 |
| GATA1 | 0.012 to 0.045 | 0.0005 to 0.073 | 0.512 to 0.540 |
| JAK1 | 0.62 to 3.26 | 0.0005 to 0.0015 | 0.780 to 0.970 |

Seven of the sixteen points carry between-variant signal above the permutation null, where
the old table had three non-zero points and no test at all.

**Spearman(rho2, ceiling) = 0.987** across the 16 points, spanning four genes, four
technologies and five depths. The ceiling rises monotonically with the estimated
signal-to-noise, which is the collapse the Methods derivation predicts and which the old
table could not test because its x-axis was constant.

`rho2` also scales with depth as the derivation requires, `eta2` being proportional to
`1/m`: for JAK1, `rho2 / m` is 0.025, 0.021 and 0.033 at m = 25, 50 and 100.

## What has not changed, and what to be careful about

**The null result is untouched.** TP53 and KRAS carry no detectable between-variant signal
at any depth: every interval covers zero and no permutation P approaches significance.
Their ceilings stay at chance. The correction makes the measurement-limited verdict for
those two datasets stronger, not weaker, because it now rests on an unbiased estimate with
a stated interval rather than on a clipped zero.

**JAK1's points rest on very few variants.** The 4m-cell requirement leaves 14, 9 and 5
variants at m = 25, 50 and 100. The intervals are correspondingly wide (rho2 = 3.26,
interval 0.11 to 5.22 at m = 100). These points anchor the high end of the collapse and
should not be read as precise values.

**The depth trend is confounded with the cohort.** Raising `m` removes variants that no
longer qualify, so the `cohort=eligible` depth series mixes depth with a changing variant
set. The `cohort=fixed` rows hold the set constant, but only 35, 49, 15 and 5 variants
survive at every depth for TP53, KRAS, GATA1 and JAK1, so those series are underpowered
and mostly non-significant. Both are in the table; neither should be quoted alone.

**The permutation test is conservative.** The observed statistic averages 20 cell splits
and the null averages 5 (`perm_splits`), so the null distribution is wider than the
observed statistic's and the reported P values are, if anything, too large.

**`null_delta2_mean` is not zero for JAK1** (1.6 to 1.7). That is expected: the null is
exact only under the hypothesis of no signal, and when signal is present, permuting mixes
it into the within-variant term. The test remains valid; the column is reported so the
behaviour is visible rather than assumed.
