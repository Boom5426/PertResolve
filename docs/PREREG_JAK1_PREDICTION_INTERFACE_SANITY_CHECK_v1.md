# Pre-registration: `JAK1_PREDICTION_INTERFACE_SANITY_CHECK_v1`

**Frozen 2026-08-03, before the check was run.** Deliberately small. This is a
pre-submission sanity check on an existing conclusion, not a new analysis
direction. It reuses the existing JAK1 data, splits, metric and stored
predictions, trains nothing, adds no dataset, and produces no new figure.

## Why

JAK1 is the manuscript's only positive example of *measurement-resolved but
model-limited* (manuscript lines 179, 202, 260; SI line 130). The TP63 screen
showed that a model sitting at chance is not automatically evidence about
features: if the estimator can only emit points inside the span of the training
variants' observed profiles, and that span is small, then "model at chance" is
forced by geometry. That failure mode has not been checked for JAK1.

**Scope statement, binding.** Model-class reachability is a *validity check
inside prediction*. It is **not** a fourth evaluation layer. The paper's spine
stays `detection -> identification -> prediction`. Nothing here may be promoted
into the main framing.

## What is already known, and does not need re-deriving

Established by inspection of the compute workspace before freezing:

- JAK1: 26 variants + WT, 5,633 cells, 2,000 genes. Splits used are
  split1/2/3/5 with train/test = 17/9, 13/13, **6/20**, 9/17, i.e. 59 scored
  (split, variant) cells. Splits are overlapping re-partitions of the same 26
  variants, not independent replicates.
- The metric is `harness.canonical_pds_cos`, tie-aware mid-rank PDS-cosine,
  candidate set = train + test of that split, chance 0.5. Canonical scoring is
  15-seed averaged with a 2,000x bootstrap over held-out variants.
- Published JAK1 numbers: oracle ceiling **0.792 (0.742-0.838)**,
  best in-house model **0.517 (Lasso-theta)**, quoted as 0.52.
- The predictor set is **mixed**. Measured as the fraction of each stored
  prediction lying outside the span of the training deltas:
  Ridge / KNN / RF / GBoost / Gene-mean / WT-null = **0.00000** (exact linear
  smoothers); Lasso-theta 0.068-0.317 and MLP-theta 0.109-0.306 (nearly);
  scGen 0.23-0.79, Biolord 0.31-0.88, CellFlow 0.29-0.34 (genuine decoders);
  PerturbNet 0.51-0.70 but confined to a 50-dim WT-PCA subspace.
- **Top-1 is never claimed for JAK1** anywhere in the manuscript or SI. Every
  JAK1 model-side claim is a PDS rank-percentile or a binary two-sample AUROC.
  The word "identifiable" applied to JAK1 is a measurement-level criterion
  (nearest sibling beyond own replicate noise), not top-1 correctness.

The concern is therefore precise rather than general: the sentence at manuscript
line 202, *"the best in-house model stayed near chance (0.52)"*, is about exactly
the family that is span-limited, and JAK1 split3 gives that family a
**6-dimensional** reachable subspace inside 2,000-dim gene space while
contributing 20 of the 59 scored cells.

## The three checks, and nothing else

Everything mirrors `oracle_ceiling.py` and `score_definitive.py` exactly: same
splits, `NSUB = 300`, `NSEED = 15`, `NBOOT = 2000`, same tie-aware PDS, candidate
set = train + test, only test variants scored. No metric is redefined.

**C1. Constructive positive control.** Feed the scorer predictions that provably
carry allele identity and confirm it says so.
- `truth`: the prediction is the variant's own observed delta. Must score
  PDS = 1.0 and top-1 = 1.0. Anything else means the scoring path is broken and
  the rest of this document is void.
- `truth + noise` at several signal-to-noise ratios: must degrade smoothly.
- `splithalf`: an independent disjoint half of the same variant's own cells,
  which is the published ceiling construction. Must reproduce ~0.792.

**C2. Model-class reachability.** For each split and each held-out variant,
compute the best prediction achievable *inside* the span of that split's training
deltas: the projection of the true delta onto that span, which maximises cosine
similarity and therefore upper-bounds every exact linear smoother
(Ridge, KNN, RF, GBoost, Gene-mean) and closely bounds Lasso and MLP. Score it
through the same pipeline. Report **per split**, because split3 is the binding
case. Independently re-derive the outside-span fraction of every stored
prediction rather than taking the pre-freeze inspection on trust.

**C3. PDS and top-1 reported separately.** Compute both for the measurement
ceiling, the span oracle, the best in-house model and the best decoder. No claim
in the manuscript depends on top-1; the purpose is to state explicitly which of
the two is reachable, so the distinction cannot be blurred later.

## Verdict rules (binding)

Let `O` be the span-oracle PDS with a bootstrap interval, and `C = 0.792` the
published ceiling.

- **`INTERFACE_OK`** if C1 passes and the span oracle's bootstrap lower bound
  exceeds 0.5 by a margin that leaves real headroom, taken as
  `O_lo > 0.5 + 0.25 * (C - 0.5) = 0.573`. The in-house family could have scored
  well above chance and did not. **Keep the model-limited conclusion**, add one
  sentence to Methods or SI.
- **`PDS_OK_TOP1_UNREACHABLE`** if the span oracle clears that bar on PDS but its
  top-1 is at or below chance. Keep the ranking wording; do not let any
  exact-identification wording enter. Since no top-1 claim currently exists, this
  outcome changes nothing but is recorded.
- **`INTERFACE_UNREACHABLE`** if `O_lo <= 0.573`, i.e. the in-house family could
  not have cleared chance by a useful margin whatever the features. Then the
  line-202 sentence about the best in-house model must be qualified, and the
  model-limited verdict rests only on the decoder arm, whose reachability is not
  in question. If the decoder arm is also at chance the overall conclusion
  survives; if the decoder arm is the only support, the sentence must say so.

A split-level version of the same rule is reported for split3 specifically, since
it is the case most likely to fail and it carries a third of the evidence.

## What this check may not do

It may not retrain anything, add a method, change a metric, change a split,
introduce a new dataset, or alter any number outside the JAK1 rows it examines.
If the outcome is `INTERFACE_UNREACHABLE` the remedy is a wording change plus a
recorded caveat, not a new experiment.

Seed 0. Runs in the `Agent` env on the remote workspace; reads
`unified/real_deltas.npz`, `allele_perturb_bench.csv`, `unified/preds5/*.npz`
and `harness.py`. Writes only new files under `unified/`.
