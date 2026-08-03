# Result: `JAK1_PREDICTION_INTERFACE_SANITY_CHECK_v1`

Run 2026-08-03 against `docs/PREREG_JAK1_PREDICTION_INTERFACE_SANITY_CHECK_v1.md`.
Code `results/reviewer_controls/jak1_interface_check.py`, outputs
`results/reviewer_controls/jak1_interface/`.

## Verdict

```
VERDICT: INTERFACE_OK
```

**The JAK1 model-limited conclusion stands.** Nothing is retracted. One Methods
paragraph was added; no number in the manuscript changed.

## The numbers

| quantity | PDS | 95% CI | top-1 |
|---|---|---|---|
| answer-carrying prediction (C1) | 1.0000 | - | 1.0000 |
| independent second measurement (published ceiling) | 0.7924 | [0.741, 0.837] | 0.247 |
| **best reachable inside the in-house family's span (C2)** | **0.6899** | **[0.638, 0.742]** | 0.049 |
| best in-house model achieved (Lasso-theta) | 0.517 | - | - |
| chance | 0.500 | - | 0.038 |

```
headroom reachable by the in-house family   65.0 %
headroom actually realised                   5.8 %
pre-registered bar (0.5 + 0.25H)            0.5730   cleared
per split, span oracle: split1 0.663  split2 0.737  split3 0.742  split5 0.617
```

C1 reproduces the published ceiling of 0.792 (0.742-0.838) to three decimals,
which is the strongest available evidence that this check scores through the same
path as `oracle_ceiling.py` and `score_definitive.py`.

## What each check settled

**C1, positive control: passed cleanly.** A prediction that carries the answer
scores PDS 1.000 and top-1 1.000; a noise ladder degrades smoothly; the split-half
arm reproduces the published ceiling. The scoring path rewards allele-specific
structure when it is present, so a null result from it is interpretable. This is
the check whose failure on TP63 invalidated an entire arm of that analysis.

**C2, model-class reachability: passed.** The in-house grid's predictions are
combinations of the training variants' measured profiles, so the projection of
each held-out variant onto that span upper-bounds the whole family. That bound is
0.690, comfortably above the 0.573 bar and above it in every split individually.
The family could have scored roughly two thirds of the way from chance to the
ceiling and scored one sixteenth of the way. **The JAK1 shortfall is about the
predictors, not about the geometry of the evaluation.**

The span restriction was verified empirically rather than assumed. Median
fraction of each stored prediction lying outside the span of the training deltas:

```
Ridge-*, KNN-*, RF-*, GBoost-*, Gene-mean, WT-null   0.00000   exact linear smoothers
Lasso-esm                                            0.023 - 0.091
MLP-esm                                              0.074 - 0.197
Lasso-theta                                          0.068 - 0.317
MLP-theta                                            0.109 - 0.306
CellFlow                                             0.294 - 0.335   genuine decoder
scGen                                                0.230 - 0.785   genuine decoder
Biolord                                              0.312 - 0.879   genuine decoder
PerturbNet                                           0.506 - 0.695   50-dim WT-PCA subspace
```

So the concern applies only to the sklearn grid, and even there the bound is not
binding. The externally trained decoders are not restricted to the span at all
and also failed, which is independent support for the same conclusion.

**C3, PDS versus top-1: reported separately, as required.** The span oracle's
top-1 is 0.0486 against a chance rate of 0.0385. The frozen rule scores that as
reachable, but the margin is negligible and the honest reading is that exact
identification is out of reach for this model class. The measurement itself only
reaches top-1 0.247. **No change is needed**, because the manuscript never claims
top-1 for JAK1: every JAK1 model-side statement is a PDS rank-percentile or a
binary two-sample AUROC, and the word "identifiable" applied to JAK1 is a
measurement-level criterion (nearest sibling beyond own replicate noise), not
top-1 correctness. That wording was already the correct one.

## One correction, logged

An intermediate run of this check reported the outside-span fraction for
Ridge / KNN / RF / GBoost as 0.17-0.65, contradicting the exact-zero expected
from the algebra. The cause was in this check, not in the benchmark: the span
basis had been built from a re-sampled per-seed delta set, whereas the stored
predictions were fitted on the fixed deterministic deltas in
`unified/real_deltas.npz`. With the correct basis the fractions are exactly
0.00000, as the algebra requires. The same error also inflated the C2 bound,
which moved from 0.654 to 0.690 once the basis was fixed. All numbers in this
document are post-correction.

## Manuscript change

One paragraph added to Methods, after the gene-shared decomposition paragraph:
`\textbf{Prediction-interface controls.}` It states the three control results and
the reachability bound, and records that top-1 is out of reach so every
model-side statement is a ranking statement. Build after the edit: rc = 0,
28 pages unchanged, 0 Overfull hbox, 0 Overfull vbox, 0 Float too large,
0 undefined, 0 LaTeX warnings.

The wording differs from the sentence proposed in planning, which said the
interface "could recover allele identity". That would overstate the result:
allele *ranking* is reachable, exact *identification* is not. The inserted text
says ranking.

## Scope, honoured

No model was retrained, no method added, no metric redefined, no split changed,
no dataset introduced. The check reused `unified/real_deltas.npz`,
`allele_perturb_bench.csv`, `unified/preds5/*.npz` and `harness.py`, and wrote
only new files. Model-class reachability remains a validity check inside
prediction; the paper's spine stays `detection -> identification -> prediction`.

## Reproduction

```bash
rsync -a alleleperturb results/reviewer_controls/jak1_interface_check.py \
      <remote>:<workspace>/_repo/...            # mirror the package layout
VCCOMPASS_BASE=<workspace> python _repo/results/reviewer_controls/jak1_interface_check.py \
      --out <workspace>/jak1_check
```

Seed 0, `NSEED = 15`, `NBOOT = 2000`, `NSUB = 300`, splits 1/2/3/5, candidate set
= train + test, tie-aware mid-rank PDS-cosine. Runs in the `Agent` env in about a
minute.
