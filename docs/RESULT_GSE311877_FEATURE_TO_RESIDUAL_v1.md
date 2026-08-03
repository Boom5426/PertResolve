# Result: `GSE311877_FEATURE_TO_RESIDUAL_PREDICTION_v1`

Run 2026-08-03 against the frozen protocol in
`docs/PREREG_GSE311877_FEATURE_TO_RESIDUAL_v1.md`. Read that file, including its
section 8 run-time addendum, before reading any number here.

## Verdict

```
VERDICT: INCONCLUSIVE_UNDERPOWERED        (both configurations)
```

TP63 **is not** established as a second model-limited system. It **is** confirmed
as a measurement-resolved system. Do not claim model-limited; do not start the
scRNA reprocessing on the strength of this.

That is the pre-registered outcome, applied mechanically. The rest of this
document is what was learned on the way, some of which is worth more than the
verdict.

---

## 1. The single most important result is a negative control, not a model score

The pre-registered pipeline was fed a feature that is predictive **by
construction** (each allele's own measured residual, reduced to 10 dimensions).
It scored near the null. That failure was then traced, and it is not a coding
error.

Every estimator in the frozen list predicts `r_hat_a = c^T Yc`: a linear
combination of the **other 17** alleles' residuals. The best any such estimator
could do is therefore bounded above by `c* = K^+ M[:, a]`, which maximises
`corr(c^T Yc, pool_a)` in closed form. That analytic bound was computed:

| configuration | residual pairwise r | **analytic model-class oracle** | measurement ceiling | reachable fraction of headroom |
|---|---|---|---|---|
| `asfrozen` | +0.824 | **PDS 0.266** | 0.861 | **-0.650** |
| `meanrm` | -0.057 | **PDS 0.806** | 0.998 | **+0.615** |

**The `asfrozen` arm is void for the model question.** Its oracle sits *below*
chance, so "every model is at chance" there is forced by geometry and carries no
information whatever about features. Had the positive control been skipped, that
arm reports every model pinned at PDS ~0.17 against a ceiling of 0.861, all
confidence intervals comfortably under the `0.25H` bar, which is a textbook
`MODEL_LIMITED_SUPPORTED` pattern and would have been completely spurious. The
`asfrozen` numbers remain valid for the measurement arm only.

Only `meanrm` can answer the question, so everything below is `meanrm`.

Second-order caveat that survives even there: the ceiling reaches **top-1 0.982**
while the model-class oracle reaches **top-1 0.019**. Rank-based PDS is partly
reachable; actual top-1 identification is structurally out of reach for this
model class at n = 18. Any future write-up must not compare a model's top-1
against the measurement's top-1 as if they were the same task.

---

## 2. Measurement ceiling

```
ceiling PDS  0.9978   95% CI [0.9935, 1.0000]   top-1 0.9815
headroom H = 0.4978    0.25H bar = 0.6245    0.50H bar = 0.7489
```

**This supersedes, and does not contradict, the 0.727 quoted in
`docs/DATASET_CANDIDATES_2026-08-03.md`.** The two use different removal bases.
The earlier audit estimated the shared-component basis from *all 18* alleles
including the one being identified, which deletes part of that allele's own
specific direction and deflates the ceiling. Here the basis is fold-local, from
the 17 training alleles only, which is the correct construction for a ceiling.
Both numbers are right under their own definition; they must never be quoted side
by side as if comparable.

---

## 3. Model results (`meanrm`, 31 feature-by-model combinations, 2,000 permutations each)

Selected rows, full table in `results/dataset_screen/predict_gse311877_main_dnp63a_meanrm.csv`.

| features | model | PDS | 95% CI | perm p | perm null | recovery vs measurement H | recovery vs oracle |
|---|---|---|---|---|---|---|---|
| position | 1-NN | 0.655 | [0.556, 0.747] | 0.069 | 0.571 | 0.311 | 0.506 |
| **esm2_window16** | **ridge** | **0.626** | [0.560, 0.697] | **0.013** | 0.536 | 0.254 | 0.413 |
| esm2_window16 | 1-NN | 0.617 | [0.517, 0.715] | 0.186 | 0.571 | 0.234 | 0.381 |
| *theta* | *mean (no features)* | *0.602* | *[0.461, 0.734]* | *1.000* | *0.602* | *0.205* | *0.333* |
| esm2_sitedelta+theta | ridge | 0.590 | [0.504, 0.676] | 0.041 | 0.520 | 0.182 | 0.295 |
| esm2_sitedelta | ridge | 0.586 | [0.501, 0.673] | 0.053 | 0.520 | 0.173 | 0.281 |
| esm2_globaldelta | ridge | 0.581 | [0.507, 0.651] | 0.091 | 0.526 | 0.162 | 0.263 |
| region | region mean | 0.560 | [0.498, 0.621] | 0.032 | 0.508 | 0.120 | 0.196 |
| theta | ridge | 0.532 | [0.484, 0.578] | 0.301 | 0.512 | 0.063 | 0.103 |

Readings that matter:

- The **top-scoring entry is `position`**, a feature that says only "where in the
  protein", not "which substitution". It beats every protein-language-model
  representation.
- The best ESM result, `esm2_window16` with ridge, is nominally significant at
  p = 0.013 against its own permutation null. **It does not survive multiple
  comparisons.** 31 combinations were run; Bonferroni needs p < 0.0016 and
  Benjamini-Hochberg at the top rank needs the same. Nothing clears either. The
  pre-registration failed to specify a correction, which is an omission recorded
  here rather than quietly ignored.
- A **feature-free baseline** that just predicts the training mean scores 0.602,
  inside the confidence interval of every model in the table. Its own residual
  correlation is 1.4e-10, i.e. it predicts nothing; its PDS is numerical
  tie-breaking noise. It is included precisely because it shows how little
  separation there is to work with at n = 18.
- `MODEL_SIGNAL_PRESENT` was not reached: it required at least two ridge feature
  sets at or above the 0.25H bar, and only one is (0.626 vs 0.6245, by 0.002).
- `MODEL_LIMITED_SUPPORTED` was not reached either: several intervals extend well
  past the bar and several permutation p values are under 0.05.

Minimum detectable PDS at 80% power with n = 18 is **0.632**, below the 0.749
half-headroom bar, so the design could in principle have seen a half-headroom
effect. But note this sits only just under the model-class oracle of 0.806, so
the usable dynamic range between "detectable" and "best possible" is thin.

---

## 4. The sharpest finding: predictions collapse onto residue position

Same-residue siblings hold position exactly constant, so only the substituted
amino acid can carry the signal. Re-run across models because the top-ranked
entry, `position/position_nn`, fails this test tautologically (with R279Q held
out, its nearest training allele *by position* is R279S).

| model | correct / 10 | mean r to own allele | mean r to sibling |
|---|---|---|---|
| position / 1-NN | 0/10 | +0.027 | +0.519 |
| **esm2_window16 / ridge** (best ESM) | **0/10** | **+0.053** | **+0.358** |
| esm2_sitedelta / ridge | 2/10 | +0.053 | +0.296 |
| esm2_sitedelta+theta / ridge | 2/10 | +0.054 | +0.304 |
| esm2_globaldelta / ridge | 0/10 | +0.026 | +0.350 |
| esm2_window16 / 1-NN | 0/10 | +0.017 | +0.401 |
| theta / ridge | 6/10 | +0.015 | -0.058 |

This is not "no signal". It is a **directional failure**: the best ESM model's
prediction is roughly seven times more similar to the *wrong* sibling than to the
true allele. The representations encode which residue was hit and not which
substitution was made, so at same-residue resolution they actively mispredict.

Statistical weight, stated conservatively: the 10 trials are 5 pairs scored
twice, so they are not independent. Treating the 5 pairs as the unit, "the model
got at least one member of the pair right" has chance 0.75 and was observed 0/5,
`p = 0.25^5 = 0.001`. The effect-size gap (+0.05 versus +0.36) is the more
informative statistic and does not depend on that counting choice. With 5 pairs
this remains a small test and should be described as such.

---

## 5. What can and cannot be said

Supported:

- TP63 GSE311877 is a **measurement-resolved** system: fold-local split-half
  identification of 18 alleles reaches PDS 0.998 and top-1 0.982 against chance
  0.500 / 0.056.
- At n = 18, **no frozen feature set predicts the allele-specific residual well
  enough to survive multiple-comparison correction**, and the single best raw
  p value comes from a local-window ESM representation, not from site-level or
  global ones.
- Protein-language-model features here **carry residue position but not
  substitution identity**, demonstrated by same-residue siblings.

Not supported, and must not be written:

- That TP63 is a second model-limited system. The pre-registered bar was not met,
  and at n = 18 it could not have been met cleanly even in principle.
- Any comparison of a model's top-1 against the measurement's top-1.
- The 0.727 and 0.998 ceilings as if they were the same quantity.

---

## 6. Consequences for the plan

Per the pre-registration's decision table, `INCONCLUSIVE_UNDERPOWERED` means:
TP63 is a measurement-ceiling system only; hold the scRNA reprocessing;
prioritise NFE2L2 access.

Final positioning of TP63, decided and frozen:

- a **measurement-resolved bulk positive control**;
- feature-to-residual prediction **underpowered**, not model-limited;
- ESM representations **encode residue position, not substitution identity**;
- **do not start the scRNA reprocessing**.

Two things changed inside that:

1. **The sibling result is worth more than the verdict.** "Representations encode
   position, not substitution" is a claim the current manuscript can make, it is
   independent of the underpowered PDS comparison, and it is exactly the
   identification-aware distinction the paper was missing. It needs a larger
   allelic series to become a headline, which is an argument *for* eventually
   getting the scRNA arm, not against.
2. **Model-class reachability belongs inside prediction, as a validity check.**
   It is **not** a fourth evaluation layer. The paper's spine stays
   `detection -> identification -> prediction`; reachability is the precondition
   that makes a prediction result readable at all, in the same way a positive
   control is. The `asfrozen` arm would otherwise have produced a clean, wrong
   `MODEL_LIMITED_SUPPORTED`, because "model at chance" is only evidence when the
   model class could have done better.

Point 2 carries one self-audit item, scoped deliberately small: the same
precondition has not been checked for the existing JAK1 model-limited claim.
That is `JAK1_PREDICTION_INTERFACE_SANITY_CHECK_v1`, a short pre-registered
check that reuses the existing JAK1 data and code, trains nothing new, and
answers only whether the original conclusion still holds. It is not a new
research direction and must not grow into one.

---

## Reproduction

```bash
# features (remote GPU, ESM2-650M)
python scripts/dataset_screen/tp63_esm2_extract.py --isoform dNp63a --out tp63_esm2_dnp63a.npz

# main run, both configurations
REMOVE_MEAN=1 python3 scripts/dataset_screen/predict_gse311877.py    # interpretable arm
REMOVE_MEAN=0 python3 scripts/dataset_screen/predict_gse311877.py    # exactly as frozen

# the controls that decide whether the above is readable at all
REMOVE_MEAN=1 python3 scripts/dataset_screen/positive_control_gse311877.py
REMOVE_MEAN=1 python3 scripts/dataset_screen/diagnose_gse311877.py
```

Seed 0, 2,000 permutations, 10,000 bootstrap resamples. Roughly 30 minutes per
configuration on 24 cores. Outputs in `results/dataset_screen/`.
