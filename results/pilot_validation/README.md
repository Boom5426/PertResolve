# 2B: prospective pilot-to-full validation

Files: `<DS>_pilot.csv` (per-perturbation pilot features + disjoint-eval rankability label),
`pilot_validation_summary.json` (LODO + calibration), `pilot_features.py`, `pilot_validate.py`.

## What was tested
The Fig 6 rankability predictor is retrospective (feature and label from the same cells). Here we
validate it PROSPECTIVELY: one permutation per perturbation sets pilot = first 50 cells, eval = next
2T cells (disjoint). Pilot effect size / pilot SNR are computed from 25-50 pilot cells; the rankability
label at depth T in {100,200} is computed by split-half S>W on the DISJOINT eval cells. A leave-one-
dataset-out logistic predictor and a training-free pilot-SNR threshold are evaluated.

## HONEST numbers (use these; adversarially verified 2026-07-10)
- **Prospective LODO AUROC, per balanced dataset (T=100): Adamson 0.85 [0.77,0.93], Norman 0.94
  [0.90,0.97], Replogle 0.99 [0.98,1.00]; mean 0.93.** Holds at 25-cell pilot (mean 0.92).
- **Negative control** (within-dataset label shuffle): null AUROC 0.49 (true 0.94 unreachable, p=0.0000).
- **Mechanistic pilot-SNR (no training), per dataset: 0.84 / 0.93 / 0.98.**
- Calibration near-diagonal (out-of-fold predicted prob -> observed rankable rate 0.04/0.00/0.22/0.79/1.00).
- Allele genes correctly at the floor: TP53 0/94, KRAS 0/96, GATA1 1/183, JAK1 7/7 rankable at T=100
  (single-class -> AUROC not evaluable; correctly predicted negatives/positives).

## DO NOT report (pooling-inflation red line, caught by adversarial verification)
- The POOLED mechanistic AUROC 0.978 is INFLATED by cross-dataset base-rate separation (positive
  fractions span 0.00-1.00); it drops to 0.71 after within-dataset rank normalization. Report PER-DATASET
  only, never the pooled value.
- allele_GATA1 LODO AUROC 0.978 is a single-positive (1/183) rank statistic, not discrimination -
  excluded from the headline.
