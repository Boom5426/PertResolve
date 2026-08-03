# Pre-registration: `GSE311877_FEATURE_TO_RESIDUAL_PREDICTION_v1`

**Frozen 2026-08-03, before any feature was computed or any model was fitted.**
Written first on purpose: the verdict rules below decide the outcome, and they
must not be reachable from the result. Any later change to this file must be
recorded as `v2` with the reason, not edited in place.

Companion: `docs/DATASET_CANDIDATES_2026-08-03.md` (why GSE311877 is a candidate),
`scripts/dataset_screen/` (the measurement-side audit this builds on).

## Question

Not "is TP63 measurable" (settled: split-half residual PDS 0.727 at chance 0.500,
permutation p = 0.0005). The question is:

> Can the reproducible, allele-specific TP63 residual be predicted from protein
> variant features?

TP63 becomes a second **model-limited** system only under
`measurement ceiling >> chance` **and** `model ~ chance`. Either of the other two
outcomes is a different, still-useful answer, and is pre-committed to below.

---

## 1. Sequence and numbering (resolved before freezing)

The 18 arrayed labels use neither UniProt isoform's numbering directly. A single
constant offset reconciles them, verified by requiring the wild-type residue of
every label to match:

```
label position N  =  TAp63-alpha (Q9H3D4, 680 aa) residue N+39
                  =  dNp63-alpha (Q9H3D4-2, 586 aa) residue N-55
18/18 wild-type residues match under both. (680 - 586 = 94 = 39 + 55.)
```

**Reference used: dNp63-alpha (Q9H3D4-2)**, because the assay is a
fibroblast-to-keratinocyte conversion and dNp63-alpha is the keratinocyte
isoform. The two sequences are identical from TAp63-alpha residue 115 onward and
every mutated site lies at 318-595, so **site-local features are
isoform-invariant here**; only whole-sequence mean embeddings could differ.
TAp63-alpha is therefore run as a declared sensitivity check, not as a separate
result.

Domain assignment from UniProt on TAp63-alpha coordinates: DNA binding 170-362,
SAM 541-607. All 8 DBD-side mutants fall in 170-362; all 10 SAM-side mutants fall
in 541-607.

---

## 2. Prediction target: fold-local residual

Raw allele delta is rejected as the primary target because a single severity
scalar already reaches PDS 0.836 on it. The target is the residual after the
shared axes are removed.

Per leave-one-allele-out fold with held-out allele `a` and training set `T`
(the other 17):

1. **Gene filter, fold-local.** Keep genes with CPM >= 1 in >= 50% of libraries,
   computed over the 4 WT libraries and the 4x17 training-allele libraries only.
   Allele `a`'s libraries do not enter the filter.
2. **Delta.** `delta_i = mean(log2CPM of the 4 reps of i) - mean(log2CPM of the 4 WT reps)`.
   WT is a control, never held out, and is used identically in every fold.
3. **PC basis, fold-local.** `U` = top **5** left singular vectors of the
   gene-space matrix of the 17 **training** deltas, column-centered. `k = 5` is
   fixed in advance from the measurement-side audit; it is not tuned.
4. **Residual.** `r_i = delta_i - U U^T delta_i`, applied to all 18 alleles using
   the training-derived `U`. The held-out allele contributes nothing to `U`.
5. **Feature standardization, fold-local.** Mean and scale from `T` only.
6. **Hyperparameters, fold-local.** Selected by an inner leave-one-out loop
   *inside* `T`. The held-out allele never influences hyperparameter choice.

Allele `a` therefore does not touch gene filtering, the residual basis,
standardization, or hyperparameter selection.

---

## 3. Features (frozen list, no additions after the fact)

Mutation-aware ESM2-650M (`esm2_t33_650M_UR50D`), final-layer per-residue hidden
states, mutant sequence built from the dNp63-alpha wild type:

| id | definition | dim |
|---|---|---|
| `esm2_sitedelta` | `h_mut[p] - h_wt[p]` at the mutated residue | 1280 |
| `esm2_window16` | mean over `p +/- 16` of `h_mut - h_wt` | 1280 |
| `esm2_globaldelta` | `mean(h_mut) - mean(h_wt)` over all residues | 1280 |
| `theta` | repo `alleleperturb.features.compute_theta` | 6 |
| `esm2_sitedelta + theta` | concatenation, each block standardized | 1286 |
| `position` | residue number | 1 |
| `region` | DBD vs SAM indicator | 1 |

Because 18 points cannot support 1280 free dimensions, every ESM block is reduced
by PCA **fitted on the training alleles only**, to `min(10, |T| - 1)` components,
inside each fold.

`theta` is `[d_hydro, d_vol, d_charge, fold_core, cat_switch, is_hotspot]`. For
TP63, `fold_core = 1` for all 18 (all sites are inside an annotated domain) and
`cat_switch = is_hotspot = 0` for all 18: no external hotspot list is applied,
because any hotspot definition drawn from phenotype data would be circular. So
`theta` carries **three** varying dimensions here. Stated, not hidden.

`||delta||` severity is recorded as a **descriptive oracle-like baseline only**.
It is computed from the held-out allele's own expression and is therefore never a
legitimate held-out input; it is reported to bound what a severity-only predictor
could reach, and is excluded from every verdict rule.

---

## 4. Models (frozen list)

Primary: **Ridge** (alpha over `logspace(-3, 6, 40)`, inner-LOO).
Secondary: **PLS** (1-3 components), **kernel ridge** (RBF).
Baselines: **feature 1-NN**, **training mean**, **region mean**,
**residue-position 1-NN**.

No MLP, transformer or GNN. At n = 18 a flexible model has no interpretive value
and would only widen the confidence interval.

---

## 5. Scoring

All predictors are scored against **one identical candidate pool**: the
**half-B residual profiles of all 18 alleles**, averaged over the three balanced
2-vs-2 replicate splits, projected with the same fold-local `U`.

- **Model**: predictor is `r_hat_a` from features.
- **Measurement ceiling**: predictor is the **half-A residual of allele `a`
  itself**, i.e. an independent measurement of the same allele.

Same pool, same metric, so the two are directly comparable.

Metrics: tie-aware **PDS** (`1 - (rank-1)/(n-1)`, midrank ties, chance 0.500),
**top-1** (chance 1/18 = 0.056), **residual Pearson** and **cosine** between
`r_hat_a` and the half-B residual of `a`, and

```
headroom recovery = (PDS_model - 0.5) / (PDS_ceiling - 0.5)
```

Uncertainty: percentile bootstrap over the 18 alleles, 10,000 resamples, for
every PDS. Null: 2,000 permutations of the allele-to-feature assignment,
re-running the full LOO loop each time.

### Secondary tests

**B. Same-residue sibling test.** R279Q/S, R304Q/T, L514D/F, C522D/G, L531E/R.
Position is held exactly constant, so only the substituted amino acid can carry
the signal. Scored as the fraction of the 10 sibling members whose predicted
residual is closer to its own half-B residual than to its sibling's. Chance 0.5;
exact binomial CI on 10 trials.

**C. Region-blocked test.** Train on DBD (8), test on SAM/TID (10), and the
reverse. **Exploratory only**, declared underpowered in advance, and it does not
enter any verdict rule.

---

## 6. Verdict rules (binding)

Let `C` = ceiling PDS point estimate, `H = C - 0.5` the headroom, and let every
model's PDS carry a 95% bootstrap interval `[lo, hi]`.

**`MODEL_SIGNAL_PRESENT`** requires all of:
- ceiling bootstrap `lo > 0.5`;
- the best primary model has `lo > 0.5` and permutation `p < 0.05`;
- its point estimate recovers `>= 0.25 H`, i.e. `PDS >= 0.5 + 0.25H`;
- at least two of the frozen feature sets agree in direction.

**`MODEL_LIMITED_SUPPORTED`** requires all of:
- ceiling bootstrap `lo > 0.5` (the measurement does resolve the allele);
- for **every** frozen model, the bootstrap **upper** bound satisfies
  `hi < 0.5 + 0.25H`, i.e. even the optimistic end cannot recover a quarter of
  the headroom;
- permutation `p > 0.05` for every model;
- a power statement: the minimum PDS effect detectable at 80% power with n = 18,
  estimated by resampling, is reported and is **below** `0.5 + 0.5H`, so the
  design could have seen a half-headroom effect had one existed.

**`INCONCLUSIVE_UNDERPOWERED`** otherwise, and specifically whenever any model's
interval covers both `0.5` and `0.5 + 0.5H`.

`p > 0.05` alone is never sufficient for `MODEL_LIMITED_SUPPORTED`. Absence of
evidence is reported as `INCONCLUSIVE_UNDERPOWERED`, which is the expected
outcome at n = 18 and is a legitimate, actionable result.

---

## 7. What each verdict triggers

| Verdict | Claim | Next action |
|---|---|---|
| `MODEL_SIGNAL_PRESENT` | TP63 is measurement-resolved **and partially predictable**; not a second model-limited system | Do **not** start scRNA reprocessing. Use TP63 as a positive control showing the framework separates measurement-limited, model-limited and predictable regimes |
| `MODEL_LIMITED_SUPPORTED` | TP63 is a second independent model-limited system | Into the main text; then commission **one** scRNA pool, not a blind full reprocess |
| `INCONCLUSIVE_UNDERPOWERED` | TP63 is a measurement-ceiling system only | Do not claim model-limited. Hold scRNA reprocessing. Prioritise NFE2L2 access |

Seed 0. numpy and pandas only for the modelling; torch plus `esm` on the remote
4090 for embedding extraction.

---

## 8. Run-time addendum, logged 2026-08-03 before any verdict was read

One mis-specification in section 2 was found while smoke-testing, and is recorded
here rather than edited away.

Section 2 step 3 takes `U` from the **allele-centred** training delta matrix. That
removes the top axes of *variation* but leaves the **grand mean response
direction** in the residual, so residuals still share a large common component
(observed pairwise r ~ 0.84). Consequence for the model arm only: any predictor
that is a convex combination of the other alleles' residuals carries that shared
offset and is actively pushed away from the held-out allele, so its PDS is **not
centred on 0.5**. The empirical null is ~0.17, confirmed by the training-mean
baseline, whose permutation null equals its observed value exactly. Section 6's
thresholds are written against `chance = 0.5` and therefore do not apply to the
model arm as literally frozen. The ceiling arm is unaffected: its predictor is an
independent measurement of the held-out allele itself.

Resolution, decided before any result was inspected: run **both** and report both.

- `asfrozen` runs section 2 exactly as written.
- `meanrm` adds the training-mean direction to the removed subspace (6 removed
  dimensions rather than 5), which is what "remove the shared axes" was meant to
  do and restores `chance = 0.5`, the assumption section 6 needs.

Neither is selected on the basis of which gives a preferable verdict. Both
verdicts are reported. Where they disagree, the disagreement itself is the
finding and the honest label is `INCONCLUSIVE_UNDERPOWERED`.

Also logged: section 5's phrase "averaged over the three balanced splits" is
ambiguous. Averaging the *profiles* would put replicates 2,3,4 in every pool
member while the ceiling's half-A predictor also holds replicates 2,3,4,
overlapping predictor with pool and inflating the ceiling. The run averages the
*scores*, the only reading that keeps half-A and half-B disjoint.
