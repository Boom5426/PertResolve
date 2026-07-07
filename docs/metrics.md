# Metric Definitions

## Ranking metrics

### PDS (Perturbation Discrimination Score)
Given a predicted delta for variant v, PDS measures whether v's measured delta
ranks highest (by distance to the prediction) among all candidate variants.
PDS = 1.0 means perfect ranking; PDS = 0.5 is chance. We report three distance
variants: cosine (PDS_cos), L1 (PDS_L1), and L2 (PDS_L2), following the
cell-eval standard (Arc Institute).

**Tie-aware scoring**: when pred_delta = zero vector, all distances are equal,
and PDS = 0.5 (not biased by list order).

## Direction metrics

### Pearson-Δ (PDCorr)
Pearson correlation between predicted and measured pseudobulk delta across all
genes. Measures global direction alignment. Used by STATE, SCALE, OCOO-T.

### Pearson-Δ̂20
Same as Pearson-Δ but restricted to the top-20 highest-variance genes across
training deltas. Introduced by SCALE to test whether signal concentrates on a
few genes.

### delta_cosine
Cosine similarity between predicted and measured delta vectors (per variant).

## DE fidelity metrics

### DE Overlap
Fraction of the top-k real DEGs (by |t-statistic|) that appear in the predicted
top-k (by |pred_delta|). Default k = 50.

### DE-LFC-Spearman
Spearman rank correlation of real LFC vs predicted delta, restricted to top-k
real DEGs. Measures effect-size ordering fidelity.

### Direction Agreement
Fraction of top-k real DEGs where the predicted sign matches the real sign.

## Reconstruction metric

### MAE
Mean absolute error of the predicted pseudobulk delta.

## Aggregation

Metrics are aggregated as **mean-of-gene-means**: per-variant scores → per-gene
mean → arithmetic mean across genes. This prevents genes with more variants
(GATA1: 255) from dominating the overall score.
