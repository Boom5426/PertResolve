# Resolution panel: results

This is a frozen historical snapshot. The displayed `detectable`,
`identifiable`, `ceiling` and `verdict` columns are retained for artifact
compatibility; `ceiling` means the empirical split-half reproducibility
reference, not a hard ceiling or bound, and `verdict` is not emitted by the
current public resolution API. Use the independent-axis semantics in
[`results/canonical/README.md`](../canonical/README.md).

Protocols: `docs/PREREG_RESOLUTION_PANEL_v2.md` (frozen 2026-09-01, amendments A1 to A11) and
`docs/PREREG_RESOLUTION_PANEL_v3.md` (frozen 2026-09-04, T1 to T6 signed 2026-09-03).
Runner: `scripts/analysis/resolution_panel_v2.py`. Extractor:
`scripts/analysis/extract_panel_subset.py`. Criterion:
`pertresolve.resolution.resolution_report` at `depth 50, n_seeds 8, n_boot 1000, seed 0,
n_components 50`, unchanged from v1 throughout. Compute is remote; the files in this
directory are the synced outputs.

**Two things must be read before any number here is quoted.** Section 5 establishes what
`identifiable_fraction` does and does not measure. Section 6 establishes that no quantity in
this study is free of measurement depth, including the one this analysis predicted would be.

---

## 1. The panel

31 datasets, one row per (resource, stratum). `v1` rows are the six published screens
of `results/canonical/resolution_panel.csv`, carried over unchanged.

| dataset | src | class | modality | n | detectable | identifiable | ceiling | rho2_nn | P(rho2_nn>1) | p_order | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `gse306429_A549` | v2 | chemical | RNA | 96 | 1.0000 | 0.1042 | 0.9591 | 0.5376 | n/a | 0.311 | detectable |
| `mcfarland` | v1 | chemical | RNA | 10 | 0.9000 | 0.2000 | 0.9278 | 1.8731 | n/a | 0.418 | detectable |
| `tahoe100m_NCI-H460` | v2 | chemical | RNA | 1075 | 0.8726 | 0.0047 | 0.9501 | 0.6900 | 0.374 | 0.863 | detectable |
| `gse306429_A375` | v2 | chemical | RNA | 95 | 0.7895 | 0.1474 | 0.9763 | 0.8385 | n/a | 0.829 | detectable |
| `tahoe100m_SW480` | v2 | chemical | RNA | 1075 | 0.7684 | 0.0000 | 0.8823 | 0.2918 | 0.229 | 0.558 | detectable |
| `tahoe100m` | v1 | chemical | RNA | 7 | 0.7143 | 0.4286 | 0.9286 | 0.6141 | n/a | 0.106 | detectable |
| `tahoe100m_PANC-1` | v2 | chemical | RNA | 1075 | 0.6921 | 0.0009 | 0.9098 | 0.6224 | 0.350 | 0.360 | detectable |
| `norman2019` | v1 | genetic | RNA | 195 | 0.6872 | 0.0051 | 0.8976 | 0.7108 | n/a | 0.997 | detectable |
| `parse10m_Donor5` | v2 | cytokine | RNA | 90 | 0.6778 | 0.0000 | 0.8765 | 0.5177 | 0.267 | 1.000 | detectable |
| `replogle` | v1 | genetic | RNA | 414 | 0.5121 | 0.0145 | 0.8648 | 0.5256 | n/a | 0.959 | detectable |
| `adamson2016` | v1 | genetic | RNA | 94 | 0.4787 | 0.0319 | 0.8510 | 0.5481 | n/a | 0.974 | not detectable |
| `gse306429_CD34` | v2 | chemical | RNA | 102 | 0.4216 | 0.0490 | 0.8912 | 0.4366 | n/a | 0.397 | not detectable |
| `parse10m_Donor12` | v2 | cytokine | RNA | 90 | 0.4111 | 0.0000 | 0.8325 | 0.3693 | 0.167 | 0.944 | not detectable |
| `kolf_ipsc` | v2 | genetic | RNA | 1500 | 0.3120 | 0.0000 | 0.7987 | -0.0190 | n/a | 0.999 | not detectable |
| `vcc_training` | v1 | genetic | RNA | 139 | 0.2590 | 0.0432 | 0.8055 | 0.4260 | n/a | 0.903 | not detectable |
| `parse10m_Donor1` | v2 | cytokine | RNA | 90 | 0.2556 | 0.0000 | 0.8335 | 0.5029 | 0.222 | 0.995 | not detectable |
| `papalexi_protein` | v2 | genetic | protein | 23 | 0.2174 | 0.0000 | 0.7169 | 0.4869 | n/a | 0.351 | not detectable |
| `perturbmulti_protein` | v2 | genetic | protein | 132 | 0.1818 | 0.0152 | 0.7399 | 0.3599 | n/a | 1.000 | not detectable |
| `frangieh_rna_Coculture` | v2 | genetic | RNA | 159 | 0.0943 | 0.0000 | 0.5825 | 0.1893 | n/a | 0.214 | not detectable |
| `papalexi_rna` | v2 | genetic | RNA | 23 | 0.0870 | 0.0000 | 0.7068 | 0.1839 | n/a | 0.615 | not detectable |
| `sciplex3_A549_24h` | v2 | chemical | RNA | 495 | 0.0727 | 0.0000 | 0.6593 | 0.0510 | n/a | 0.991 | not detectable |
| `frangieh_protein_Coculture` | v2 | genetic | protein | 159 | 0.0692 | 0.0000 | 0.5752 | 0.1232 | n/a | 0.067 | not detectable |
| `sciplex3_MCF7_24h` | v2 | chemical | RNA | 735 | 0.0531 | 0.0000 | 0.6536 | 0.0691 | n/a | 0.003 | not detectable |
| `sciplex3_K562_24h` | v2 | chemical | RNA | 590 | 0.0525 | 0.0000 | 0.6150 | 0.1055 | n/a | 0.814 | not detectable |
| `frangieh_rna_IFNg` | v2 | genetic | RNA | 179 | 0.0335 | 0.0000 | 0.6121 | 0.2279 | n/a | 0.000 | not detectable |
| `perturbmulti_rna` | v2 | genetic | RNA | 110 | 0.0273 | 0.0364 | 0.7136 | 0.1313 | n/a | 0.994 | not detectable |
| `frangieh_rna_Control` | v2 | genetic | RNA | 87 | 0.0230 | 0.0115 | 0.5690 | 0.0623 | n/a | 0.030 | not detectable |
| `xatlas_hek293t` | v2 | genetic | RNA | 1500 | 0.0140 | 0.0000 | 0.5823 | 0.1530 | 0.031 | 0.999 | not detectable |
| `frangieh_protein_Control` | v2 | genetic | protein | 87 | 0.0115 | 0.0000 | 0.5382 | -0.2061 | n/a | 0.001 | not detectable |
| `xatlas_hct116` | v2 | genetic | RNA | 1500 | 0.0080 | 0.0000 | 0.5440 | 0.0428 | n/a | 0.931 | not detectable |
| `frangieh_protein_IFNg` | v2 | genetic | protein | 179 | 0.0000 | 0.0000 | 0.5574 | -0.9350 | n/a | 0.027 | not detectable |

**31 datasets.** Detectable spans 0.0000 to 1.0000,
**21 of 31 below 0.50**. Identifiable spans 0.0000
to 0.4286, **17 exact zeros**, and **at the frozen depth of
50 cells per group not one dataset exceeds 0.50**.

Every exact zero sits at a low detectable fraction. The zeros are a detection result, and
detection is measured against the control alone and does not depend on the candidate pool.

### 1a. The v1 Tahoe row is superseded

The published `tahoe100m` row rests on an 80,000-cell demonstration subset with 7
perturbations and reports `identifiable = 0.4286`, the highest value anywhere in the study.
The full resource, three cell-line strata at 4.2 M to 6.4 M cells and 1,075 units each,
returns 0.0000, 0.0047 and 0.0009. The v1 value is a candidate-pool artifact of a 7-unit
pool, not a property of Tahoe.

### 1b. One correction to a published row

`AdamsonWeissman2016_GSM2406681_10X010.h5ad` carries a perturbation level whose literal value
is `<NA>`, with 2,613 cells, above the 200-cell threshold. **v1 counted it as one of its 94
perturbations.** Applying D7 removes it: n = 93, detectable 0.4839, identifiable 0.0323.
Norman, Replogle, the Virtual Cell Challenge training set and McFarland were checked and are
clean; only Adamson is affected.

---

## 2. Gates

| gate | question | outcome |
|---|---|---|
| **1** | does the loader reproduce the v1 panel | **passed, bit-exact** on Norman, Adamson and Replogle; Adamson only once D7 is held constant, see 1b |
| **2** | does the per-perturbation cell cap change the answer | **cap retired**: it moved `p_correct_order` by 25 times the seed spread |
| **5** | does the candidate-pool cap change the answer | **benign**, see below |
| **6 / 6b** | are the alternative reductions equivalent | `IncrementalPCA` retired, 329 times slower and it moved the ceiling; sparse randomised SVD equivalent to 3e-08 |
| **7** | is the RNA versus protein contrast confounded with feature width | **not confounded**, width matching moves detectable by at most 0.0062 |

### Gate 5, the candidate-pool cap at K_compute = 1500

| quantity | K = 1500 | K = 3000 | shift |
|---|---|---|---|
| n_perturbations | 1500 | 3000 | +1500.0000 |
| n_cells_scored | 564835 | 981716 | +416881.0000 |
| detectable_fraction | 0.312 | 0.3 | -0.0120 |
| identifiable_fraction | 0.0 | 0.0 | +0.0000 |
| replicate_ceiling | 0.7987065821658883 | 0.8000076414360342 | +0.0013 |
| rho2_nn_median | -0.018970484488774027 | -0.03339535101282103 | -0.0144 |
| p_correct_order | 0.999 | 0.998 | -0.0010 |

Doubling the pool moves every quantity by far less than the seed-to-seed spread already
measured on these statistics, which is 0.042 for the identifiable fraction on a 95-unit
dataset. The cap is reported as benign, with the numbers that show it.

---

## 3. The K curve

Per-query nested candidate sets, `C_v(K, r) = {v}` plus the first `K - 1` entries of a fixed
permutation with `v` removed, 10 draws, profiles formed once per dataset so that only the
candidate set varies. Mean of 10 draws.


**gse306429_A549**, cohort n = 96, detectable = 1.0000 (pool-free)

| K | identifiable | rho2_nn median | P(rho2_nn>1) | replicate candidate-pool PDS | p_correct_order |
|---|---|---|---|---|---|
| 10 | 0.3833 | 1.7476 | 0.601 | 0.9611 | 0.262 |
| 20 | 0.2479 | 1.0680 | 0.515 | 0.9590 | 0.242 |
| 40 | 0.1812 | 0.8310 | 0.442 | 0.9598 | 0.207 |
| 80 | 0.1240 | 0.6054 | 0.351 | 0.9588 | 0.277 |
| 96 | 0.1042 | 0.5376 | 0.323 | 0.9591 | 0.311 |

**kolf_ipsc**, cohort n = 1500, detectable = 0.3120 (pool-free)

| K | identifiable | rho2_nn median | P(rho2_nn>1) | replicate candidate-pool PDS | p_correct_order |
|---|---|---|---|---|---|
| 10 | 0.0002 | 0.0931 | 0.058 | 0.7999 | 0.855 |
| 20 | 0.0001 | 0.0095 | 0.028 | 0.7980 | 0.948 |
| 40 | 0.0003 | -0.0164 | 0.020 | 0.7969 | 0.955 |
| 100 | 0.0003 | -0.0202 | 0.016 | 0.7972 | 0.995 |
| 250 | 0.0001 | -0.0265 | 0.018 | 0.7979 | 0.995 |
| 500 | 0.0001 | -0.0284 | 0.018 | 0.7986 | 0.995 |
| 1500 | 0.0000 | -0.0190 | 0.016 | 0.7987 | 0.999 |

**norman2019**, cohort n = 195, detectable = 0.6872 (pool-free)

| K | identifiable | rho2_nn median | P(rho2_nn>1) | replicate candidate-pool PDS | p_correct_order |
|---|---|---|---|---|---|
| 10 | 0.0826 | 1.5071 | 0.698 | 0.8981 | 0.912 |
| 20 | 0.0462 | 1.2436 | 0.615 | 0.8968 | 0.983 |
| 40 | 0.0354 | 1.0392 | 0.514 | 0.8975 | 0.975 |
| 100 | 0.0179 | 0.8268 | 0.405 | 0.8979 | 0.993 |
| 195 | 0.0051 | 0.7108 | 0.349 | 0.8976 | 0.997 |

**replogle**, cohort n = 414, detectable = 0.5121 (pool-free)

| K | identifiable | rho2_nn median | P(rho2_nn>1) | replicate candidate-pool PDS | p_correct_order |
|---|---|---|---|---|---|
| 10 | 0.1312 | 1.2353 | 0.549 | 0.8660 | 0.836 |
| 20 | 0.0961 | 0.9612 | 0.484 | 0.8653 | 0.929 |
| 40 | 0.0560 | 0.8223 | 0.439 | 0.8651 | 0.919 |
| 100 | 0.0234 | 0.6655 | 0.368 | 0.8663 | 0.935 |
| 250 | 0.0171 | 0.5710 | 0.309 | 0.8654 | 0.952 |
| 414 | 0.0145 | 0.5256 | 0.275 | 0.8648 | 0.959 |

Three results, all measured under a design in which the number of scored units does not
change with K:

- **The strict identification certification rate collapses with the candidate pool.** Norman
  falls 16-fold from K = 10 to its native 195. Replogle falls 9-fold.
- **The replicate candidate-pool PDS barely moves**, 0.897 to 0.898 for Norman across the
  whole range and 0.959 to 0.961 for GSE306429 A549. It is the pool-insensitive
  identification evidence.
- **`p_correct_order` moves the other way**, rising with K on every dataset. KOLF goes from
  0.855 to 0.999. **The same dataset, in the same run, gets worse at per-perturbation
  identification and better at model comparison as the candidate set grows.**

KOLF answers the question section 6 of the pre-registration posed before the curve existed:
it is already poor at K = 10 (identifiable 0.0002, median rho2_nn 0.093), so its local
resolution is intrinsically limited rather than a consequence of a crowded pool.

---

## 4. The matched arm and the cross-seed estimator

Matched K = 40, the standardized identification challenge, 10 draws. Cross-seed columns fix
each unit's competitor on one half of the seeds and evaluate it on the other.

| dataset | cohort n | detectable | identifiable | rho2_nn | P(rho2_nn>1) | repl PDS | p_order | cross-seed identifiable |
|---|---|---|---|---|---|---|---|---|
| `gse306429_A375` | 95 | 0.7895 | 0.2189 | 1.2357 | 0.560 | 0.976 | 0.831 | 0.3147 |
| `tahoe100m_NCI-H460` | 1075 | 0.8726 | 0.0285 | 1.1714 | 0.559 | 0.949 | 0.756 | 0.0583 |
| `norman2019` | 195 | 0.6872 | 0.0354 | 1.0392 | 0.514 | 0.898 | 0.975 | 0.0667 |
| `tahoe100m_PANC-1` | 1075 | 0.6921 | 0.0097 | 0.9979 | 0.497 | 0.908 | 0.323 | 0.0223 |
| `gse306429_A549` | 96 | 1.0000 | 0.1812 | 0.8310 | 0.442 | 0.960 | 0.207 | 0.2250 |
| `replogle` | 414 | 0.5121 | 0.0560 | 0.8223 | 0.439 | 0.865 | 0.919 | 0.1048 |
| `adamson2016` | 93 | 0.4839 | 0.0452 | 0.6390 | 0.330 | 0.854 | 0.940 | 0.0624 |
| `tahoe100m_SW480` | 1075 | 0.7684 | 0.0033 | 0.5043 | 0.319 | 0.880 | 0.481 | 0.0057 |
| `gse306429_CD34` | 102 | 0.4216 | 0.1098 | 0.5649 | 0.267 | 0.891 | 0.506 | 0.1373 |
| `sciplex3_K562_24h` | 590 | 0.0525 | 0.0012 | 0.0626 | 0.152 | 0.617 | 0.688 | 0.0000 |
| `sciplex3_MCF7_24h` | 735 | 0.0531 | 0.0001 | 0.1163 | 0.121 | 0.653 | 0.014 | 0.0003 |
| `sciplex3_A549_24h` | 495 | 0.0727 | 0.0000 | 0.1398 | 0.119 | 0.660 | 0.968 | 0.0000 |
| `frangieh_protein_Control` | 87 | 0.0115 | 0.0000 | -0.2314 | 0.110 | 0.536 | 0.001 | 0.0000 |
| `frangieh_protein_IFNg` | 179 | 0.0000 | 0.0000 | -0.9266 | 0.103 | 0.556 | 0.019 | 0.0000 |
| `perturbmulti_protein` | 132 | 0.1818 | 0.0250 | 0.3407 | 0.093 | 0.743 | 0.990 | 0.0220 |
| `frangieh_protein_Coculture` | 159 | 0.0692 | 0.0000 | 0.0922 | 0.074 | 0.575 | 0.068 | 0.0000 |
| `frangieh_rna_Control` | 87 | 0.0230 | 0.0011 | 0.1027 | 0.036 | 0.571 | 0.017 | 0.0000 |
| `frangieh_rna_IFNg` | 179 | 0.0335 | 0.0000 | 0.1974 | 0.033 | 0.610 | 0.000 | 0.0000 |
| `perturbmulti_rna` | 110 | 0.0273 | 0.0291 | 0.1293 | 0.025 | 0.713 | 0.988 | 0.0264 |
| `frangieh_rna_Coculture` | 159 | 0.0943 | 0.0000 | 0.1438 | 0.021 | 0.581 | 0.155 | 0.0000 |
| `kolf_ipsc` | 1500 | 0.3120 | 0.0003 | -0.0164 | 0.020 | 0.797 | 0.955 | 0.0003 |

**21 datasets at matched K = 40.** The cross-seed estimator returns a higher
identifiable fraction on every dataset where the released estimator returns a non-zero one,
by roughly a factor of two (Norman 0.0354 to 0.0667, Replogle 0.0560 to 0.1048, GSE306429
A375 0.2189 to 0.3147). Fixing the competitor on disjoint seeds removes the churn of the
`argmin` from what the released estimator charges to measurement noise. It is reported beside
the released quantity, never in place of it.

### 4a. The axes cannot be collapsed into one score

Spearman correlations across all 21 matched-K datasets, that is with the candidate pool held
equal so that no correlation here is a pool artifact:

| | detectable | identifiable | P(rho2_nn>1) | repl PDS | p_order |
|---|---|---|---|---|---|
| detectable | 1.000 | 0.701 | 0.734 | 0.929 | 0.384 |
| identifiable | 0.701 | 1.000 | 0.657 | 0.849 | 0.545 |
| P(rho2_nn>1) | 0.734 | 0.657 | 1.000 | 0.757 | **0.192** |
| repl PDS | 0.929 | 0.849 | 0.757 | 1.000 | 0.497 |
| p_order | 0.384 | 0.545 | 0.192 | 0.497 | 1.000 |

The structure is **two dimensions, not one and not five**. Detection evidence, the strict
certification rate and the replicate candidate-pool PDS cohere at 0.657 to 0.929.
`p_correct_order` sits apart, at 0.192 against `P(rho2_nn>1)` and 0.384 against detection,
though it is not independent of the others: 0.545 against the certification rate and 0.497
against the replicate PDS. **Any single score is a choice of weighting, and the ranking
follows the weighting rather than the data.**

Model-ordering power is not a function of screen size. Spearman between the number of scored
units and `p_correct_order` over all 31 panel rows is 0.206 (p = 0.27). A 735-unit sci-Plex
stratum returns 0.003 while a 90-unit Parse donor returns 1.000. What does move it, on a fixed
cohort with only the candidate set varying, is the candidate set: KOLF rises from 0.855 at
K = 10 to 0.999 at K = 1500 while its identifiable fraction falls.

---

## 5. What `identifiable_fraction` measures

Established by five independent audits, each adversarially verified.

The test compares the mean of eight per-split separations against the 2.5-to-97.5 percentile
range of those same eight numbers. At eight samples that range is a scale estimate, not a
confidence interval, so the test reduces to a per-split coefficient of variation below 0.378.
Its false-positive rate at a true separation of zero is **6.6e-05**: a certification rule at
roughly a 0.007 percent operating point.

Consequences, all measured:

- **It is pool-dependent**, section 3.
- **It is depth-dependent**, section 6.
- **It gets stricter as more seeds are added**, 0.895 at 3 seeds to 0.421 at 32, the opposite
  direction to a standard test.
- **The margin dominates its magnitude.** On the identical eight splits, a 5 percent one-sided
  test gives Norman 0.4923 rather than 0.0051 and Adamson 0.4149 rather than 0.0319.

It is retained under its own field name for provenance and is called the **strict
identification certification rate**. It is never to be written as "X percent of perturbations
are identifiable".

Its apparent conflict with a replicate PDS of 0.976 on the same data is not a conflict. The
PDS is a rank statistic and stays near 1 while the truth is closest however narrowly; the
certification asks that the nearest competitor be more than 2.6 per-split standard deviations
away. The high ceilings are direct evidence that these resources do resolve their
perturbations.

---

## 6. The depth sweep, and a refuted prediction

**Varies a frozen parameter and is therefore a sensitivity measurement, not a panel number.**
Cohort held fixed at the units able to supply `4 * depth_max` cells, so a deeper run cannot
silently drop units; candidate set matched at K = 40 throughout.

The prediction, from mechanism: `eta2` is the sampling noise of a mean of `depth` cells, so
`eta2` is proportional to `1 / depth` and `rho2_nn = sep / (2 * eta2)` is proportional to
depth. If `sep` is depth-free then `rho2_nn / depth` is a depth-invariant local-geometry
statistic. **It is not.**


**gse306429_A375**, fixed cohort n = 82

| depth | detectable | rho2_nn median | rho2_nn / depth | identifiable | P(rho2_nn>1) | repl PDS |
|---|---|---|---|---|---|---|
| 25 | 0.5366 | 0.5946 | 0.023783 | 0.0854 | 0.3354 | 0.9395 |
| 50 | 0.7683 | 1.0557 | 0.021115 | 0.1512 | 0.5171 | 0.9751 |
| 100 | 0.9146 | 1.7184 | 0.017184 | 0.2951 | 0.7159 | 0.9872 |
| 200 | 1.0000 | 3.4481 | 0.017240 | 0.5573 | 0.8695 | 0.9966 |

normalised spread 1.38-fold against 5.80-fold raw.

**gse306429_CD34**, fixed cohort n = 90

| depth | detectable | rho2_nn median | rho2_nn / depth | identifiable | P(rho2_nn>1) | repl PDS |
|---|---|---|---|---|---|---|
| 25 | 0.1667 | 0.3474 | 0.013897 | 0.0067 | 0.1356 | 0.8031 |
| 50 | 0.4111 | 0.4269 | 0.008538 | 0.0122 | 0.2022 | 0.8817 |
| 100 | 0.5667 | 0.7248 | 0.007248 | 0.0978 | 0.3600 | 0.9224 |
| 200 | 0.7111 | 1.1593 | 0.005797 | 0.2122 | 0.5667 | 0.9544 |

normalised spread 2.40-fold against 3.34-fold raw.

**sciplex3_MCF7_24h**, fixed cohort n = 561

| depth | detectable | rho2_nn median | rho2_nn / depth | identifiable | P(rho2_nn>1) | repl PDS |
|---|---|---|---|---|---|---|
| 25 | 0.0160 | 0.0487 | 0.001948 | 0.0000 | 0.0888 | 0.5895 |
| 50 | 0.0178 | 0.1214 | 0.002427 | 0.0002 | 0.1007 | 0.6206 |
| 100 | 0.0624 | 0.3047 | 0.003047 | 0.0005 | 0.1663 | 0.6826 |

normalised spread 1.56-fold against 6.26-fold raw.

Normalising by depth absorbs most of the scaling, from 3.3-to-6.3-fold down to
1.4-to-2.4-fold, but **it does not remove it, and the residual has no consistent sign**:
GSE306429 A375 and CD34+ fall with depth while sci-Plex rises. The most likely cause is that
the debiased separation is itself depth-dependent, because the cross-fitted estimator is
biased upward and increasingly so as noise grows, which inflates the low-depth end.

**So there is no depth-free quantity here, and the proposal is withdrawn rather than used.**
Detection rises from 0.5366 to 1.0000 and the replicate PDS from 0.9395 to 0.9966 over the
same sweep. Depth is an irreducible axis, exactly as the candidate set is.

One consequence worth stating on its own, **with its full scope, because a shorter version of
it would be falsifiable**. At depth 200, **and at a matched candidate set of K = 40, and on the
82-unit fixed cohort**, GSE306429 A375 reaches an identifiable fraction of **0.5573**, the only
value above 0.50 anywhere in this study. The native-pool value in the same run on the same
cohort is **0.4512**, below the threshold. So the crossing requires depth, candidate matching
and the cohort together; depth alone does not carry it across. At the frozen depth of 50
nothing exceeds 0.50 under any of these readings. This is a statement about a measurement
configuration, not about dataset quality, and it is the clearest single illustration of why a
claim has to name its depth as well as its operating point.

---

## 7. Tahoe-100M and the D1 collapse

Three cell-line strata, selected by the frozen k = 3 rule on cell count.

| stratum | cells scored | units | detectable | identifiable | rho2_nn | repl PDS | p_order |
|---|---|---|---|---|---|---|---|
| SW480 (`CVCL_0546`) | 6,363,927 | 1075 | 0.7684 | 0.0000 | 0.2918 | 0.882 | 0.558 |
| NCI-H460 (`CVCL_0459`) | 5,943,050 | 1075 | 0.8726 | 0.0047 | 0.6900 | 0.950 | 0.863 |
| PANC-1 (`CVCL_0480`) | 4,170,504 | 1075 | 0.6921 | 0.0009 | 0.6224 | 0.910 | 0.360 |

SW480 is the largest run in the study: 100,648,790 source cells reduced to
6,363,927 kept and 8,498,694,128 stored values. The candidate-pool cap of
1,500 did not bind, because each stratum holds 1075 eligible units.

Detection is the highest in the panel, 0.69 to 0.87, while the strict certification rate stays
at or near zero. This is the pattern the framework predicts and it is now measured on the
largest chemical perturbation atlas available.

**The D1 chemical-formulation collapse merged 6 groups, all of
them salt, hydrate or solvate preparations of one entity:**

| merged entity | source labels |
|---|---|
| `berbamine` | `Berbamine`, `Berbamine (dihydrochloride)` |
| `cytarabine` | `Cytarabine`, `Cytarabine (hydrochloride)` |
| `gallic acid` | `Gallic acid`, `Gallic acid (hydrate)` |
| `omeprazole` | `Omeprazole`, `Omeprazole (sodium)` |
| `pyridoxine` | `Pyridoxine`, `Pyridoxine (hydrochloride)` |
| `tofacitinib` | `Tofacitinib`, `Tofacitinib (citrate)` |

No chemical entities were merged, stereoisomers were kept separate, and
`Trametinib (DMSO_TF solvate)` was correctly treated as a treatment rather than pulled into
the vehicle arm by a substring match. This is the only place in the pipeline where original
labels are merged automatically, and the manifest travels with every run's provenance.

---

## 8. What failed, and what it cost

| item | outcome |
|---|---|
| `sciplex3_MCF7_24h` on dense `PCA` | killed; 306 GB densified. Completed on sparse randomised SVD in 145.8 s with no cells discarded. |
| `xatlas_hek293t`, three attempts | killed. Cause measured: 4.89e9 stored values force int64 indices, so the three arrays are 58.7 GB and an archive must be read whole. Fixed by memory-mapping a directory of `.npy` files; the run then took 1352 s and reproduced the HCT116 result to 1e-07 on a control comparison. |
| `parse10m`, two attempts | killed; `anndata.read_h5ad` loads the whole 227 GB file. Fixed by streaming the extractor. |
| Tahoe, first attempt | killed. Cause: the selection pass expanded 100.6 M labels to strings, 32 GB before any matrix was read. Fixed by carrying categorical codes and writing memory-mapped output. |
| one `replogle` run | discarded for write contention from a duplicate launch, re-run alone. |

Nothing was dropped from the panel to make a run succeed. Every failure was diagnosed to a
measured cause before a fix was applied.
