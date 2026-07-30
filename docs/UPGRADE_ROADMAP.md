# AllelePerturb Upgrade Roadmap (long-term, to structurally secure Nature Methods)

Author-facing planning doc. Text and tables only; figures are drawn separately by the
author from the textual figure specifications in Section 8. Companion to
`docs/REMOTE_INVENTORY.md` (what exists on the server) and the two `CLAUDE.md` files
(red lines). Status tags: DONE / IN-MS (already in manuscript) / TODO / BLOCKED.

Numbers already established are quoted with their source; every value to be produced by a
new run is written `TBD(run:<script>)` and must never be filled with a guessed number.

---

## 1. Positioning: what this paper becomes

Current center of gravity: "a protein-coding variant benchmark on which every model sits at
chance, explained by a measurement-resolution floor." That is a clean negative result but
reads to an editor as one step past Ahlmann-Eltze & Huber 2025.

Target center of gravity (the thesis the whole upgrade serves):

> Whether single-cell perturbation models can be *fairly compared at all* is set by whether
> the measurement provides a large enough perturbation-specific signal window; we define,
> derive and validate a measurement-aware criterion that predicts when a benchmark can
> recover the true ordering of models, and demonstrate it at allele resolution and on
> gene-level Perturb-seq atlases.

The upgrade is what converts an *observation* (floor exists) into a *validated criterion*
(the floor predicts benchmark validity). The negative result stays the honest backbone; we
do not chase a model that "wins."

### Hard guards (do not let the framework ambition break these)

- The null stays honest: no reframing, threshold-tuning or split selection to make any model
  beat chance. Variant-level split unit; test variants never touch fitting/normalization.
- Report per-gene before pooled. No evaluation-definition change to inflate a method.
- Every new claim traces to a committed result file.
- The new experiments are *metrology validation*, not model-performance claims. They must be
  constructed so they can only strengthen or bound the null, never launder it into a positive.
- Do not gate the paper on data we cannot get (RUNX1 reprocessing) or on a title change.

---

## 2. Prerequisites (must clear before quoting any new number)

| # | Prerequisite | Why | Status |
|---|---|---|---|
| P1 | Lock ONE canonical scorer (`unified/harness.py` tie-aware `PDS_cos`); reconcile `definitive_summary.csv` vs `unified_summary5.csv` (PerturbNet 0.542 vs 0.561; scVIDR NaN handling) | Every new (window, PDS, ranking-recovery) point must be on one metric or the collapse plots are polluted | TODO |
| P2 | Confirm `floor_audit/*_rankability.csv` S and W map exactly to the manuscript D_null / D_self definitions | The new law and tables reuse S, W; a definitional drift breaks traceability | TODO |
| P3 | Deprecate leaky artifacts: `results/second_probe_predictor_results.csv`, `fig_config.load_predictor` | Open follow-up from the honest-predictor rebuild; stale leaky CSV must not be cited | TODO |
| P4 | Fix one Biolord number (manuscript 0.49 vs inventory 0.52) once P1 is locked | Internal consistency | TODO |

---

## 3. Master plan (overview)

Tier A = required for "structurally NM". Tier B = raises reviewer-survival. Tier C = reach.
Cost is incremental compute on the 4090 reusing existing pseudobulks + harness.

| ID | Tier | Upgrade | New compute | Feeds | Depends | Status |
|----|------|---------|-------------|-------|---------|--------|
| 0  | A | Honesty/precision text convergence (Section 4) | none | abstract, Discussion | P1 not needed | TODO |
| 1A | A | Oracle measurement ceiling `PDS_oracle` | low | Fig 3 | P1 | TODO |
| 1B | A | Pairwise allele-allele discriminability | low | Fig 3 | P1,P2 | TODO |
| 1C | A | Controlled-predictor benchmark-resolution curve | medium | Fig 4 (keystone) | P1,1A | TODO |
| 1D | A | Analytical scaling law g(rho) + collapse | low | Fig 4 | P2,1C | TODO |
| 2A | B | Metric-family generalization (energy/MMD/Wasserstein/classifier-2-sample) | medium | Fig 4 / Table C | P1 | TODO |
| 2B | B | Pilot-to-full prospective validation | medium | Fig 5 | 1A | TODO |
| 2C | B | Replicate-aware noise check | low | Fig 3 note / SI | data check | BLOCKED? |
| 2D | B | Continuous calibrated reliability score | low | Fig 4 / Methods | 1C | TODO |
| 3A | C | One independent variant dataset (break gene-tech-celltype confound) | high | ext. validation | author matrices | BLOCKED |
| 3B | C | Deeper gene-level atlas panel (7-10 datasets) | medium | Fig 4 extension | scPerturb pull | TODO |
| 3C | C | Title / positioning upgrade | none | title, abstract | 1C,1D done | DONE 2026-07-30 |

Explicitly de-scoped now: adding more models (marginal value ~0; the paper is about *when a
leaderboard is valid*, not who leads it).

---

## 4. Phase 0: immediate honesty/precision text (no compute)

These are pure precision fixes, fully inside the red lines; they should ship regardless of the
science phases. Written as before/after so they can be applied verbatim.

1. Abstract causal softening.
   - Before: "the discrimination floor therefore bounds every model class we tested, not only
     simple baselines."
   - After: "none of the evaluated model classes exceeded the empirical discrimination floor
     under the available measurement regimes."
   - Reason: the original asserts the floor is the *proven cause* of failure; we have shown
     co-occurrence, not (yet) causation. 1C supplies the causation later.

2. Replogle qualifier (every occurrence).
   - Before: "roughly half of Replogle perturbations are un-rankable."
   - After: "roughly half of Replogle perturbations are un-rankable under the specified
     split-half energy-distance criterion at native depth."
   - Reason: prevent the reading "biologically undetectable."

3. Guidance scoping.
   - Before: "provides effect-size-conditioned guidance."
   - After: "provides effect-size-conditioned triage and diagnostic guidance."
   - Reason: no prospective depth calibration exists yet (2B would add it).

4. JAK1 confound caveat (Results Fig 3 + Discussion limitation).
   - Add one sentence: "Because JAK1 differs from the other genes in assay technology, cell
     line and stimulation, its wider window reflects a more favorable effect-size-to-noise
     regime for that dataset, not necessarily stronger allele biology; gene identity and
     measurement regime are confounded across the four datasets."
   - Reason: the four datasets fully confound gene x technology x cell type; the manuscript
     must not let a gene name read as a biological explanation.

5. External-model central sentence (keep the honest framing already close to this).
   - Preferred: "Across all evaluated model classes and public interfaces, no method overcame
     the empirical discrimination floor," rather than "no current model can perform
     allele-resolution prediction."

---

## 5. Phase 1: the keystone closed loop (measurement window -> benchmark reliability)

This is the scientific core. Four interlocking pieces produce one new main figure (Fig 4)
plus two panels folded into Fig 3. All reuse `grid_cbv.npz` per-variant cells and the
canonical `PDS_cos`; no new data.

### 1A. Oracle measurement ceiling `PDS_oracle`

Purpose: an achievable upper bound on PDS given finite-sample ground truth; the single
cleanest separation of *measurement limitation* from *model limitation*.

Construction:
- For gene g, held-out variant set T. For each variant v with >= 2m cells, split its cells
  deterministically into halves A_v, B_v of size m.
- Treat pseudobulk(A_v) as an "oracle prediction" of v and pseudobulk(B_v) as evaluation
  truth. Score the whole set with the canonical tie-aware `PDS_cos`, distances computed
  between {A_v} predictions and {B_v} truths.
- Average over K deterministic seeds (match the 15-subsample convention). Report per gene and
  as a function of depth m over the existing 50/100/150/300 grid.

Interpretation: any real model's PDS <= `PDS_oracle`. If `PDS_oracle` ~ 0.5 at full depth, no
model can do better and the floor is the dominant cause. If `PDS_oracle` is high while the best
model sits at 0.5, that gene has a genuine computational gap.

Inputs: `grid_cbv.npz`, `harness.PDS_cos`. Expected (hypothesis, to verify):
TP53/KRAS `PDS_oracle` near 0.5 even at native depth; JAK1 high; GATA1 intermediate.
Output file: `TBD(run:oracle_ceiling.py -> results/oracle_ceiling.csv)`.

### 1B. Pairwise allele-allele discriminability

Purpose: close the conceptual gap that `D_null` measures v-vs-WT (detection) while PDS needs
v-vs-other-alleles (identification). A variant can be easy vs WT yet unresolvable vs its
siblings.

Construction:
- Compute the full pairwise pseudobulk distance matrix `D(v_i, v_j)` per gene in the PDS
  space, and per-variant split-half noise `D_self(v)`.
- Pairwise separability: `sep(v_i, v_j) = D(v_i, v_j) / sqrt(D_self(v_i) * D_self(v_j))`.
- A pair is "resolvable" when `sep` exceeds the same style of criterion used for rankability
  (bootstrap/permutation consistent with `S > W`). Report per-gene: fraction of allele-pairs
  resolvable, distribution of `sep`, and the fraction of variants detectable-vs-WT for
  contrast.

Interpretation: expect detectable-vs-WT >> pairwise-resolvable for TP53/KRAS (variants share
gene-level direction, pairwise signal tiny), both high for JAK1. This aligns the measurement
analysis with the PDS task and explains PDS ~ 0.5 mechanistically.

Inputs: pseudobulks, `D_self` from `floor_audit`. Output:
`TBD(run:pairwise_resolvability.py -> results/pairwise_resolvability.csv)`.

### 1C. Controlled-predictor benchmark-resolution curve (KEYSTONE)

Purpose: prove the diagnostic predicts *benchmark validity*, i.e. whether the observed scores
recover the true ordering of models. This is the experiment that converts the paper from
retrospective diagnosis to validated criterion.

Construction:
- Three-way per-variant split: build-half (estimate true effect `delta_v`), truth-half
  (evaluation target), and independent noise. Using a separate build-half avoids trivial
  self-matching, so the alpha=1 predictor's ceiling is exactly `PDS_oracle` from 1A, not 1.0
  (this ties 1A and 1C together and keeps the ceiling honest).
- Family of synthetic predictors with KNOWN quality order, alpha in {0, 0.25, 0.5, 0.75, 1.0}:

      delta_hat(v, alpha) = delta_bar_g + alpha * (delta_v - delta_bar_g) + eps,   eps ~ N(0, beta^2 * Sigma_hat)

  alpha=1 -> fully allele-specific; alpha=0 -> gene-mean (direction-only); beta is an optional
  second noise knob. True quality is monotone in alpha by construction.
- Score all alpha-predictors with the canonical PDS. Measure whether observed PDS recovers the
  true alpha order:
  - Kendall tau / Spearman between alpha and observed PDS;
  - false-inversion rate (fraction of alpha_i > alpha_j pairs with PDS_i < PDS_j);
  - minimum detectable alpha-gap; bootstrap probability that the observed winner is alpha=1.
- Plot each recovery metric against the gene's measurement window (`D_self/D_null`, or rho
  from 1D, or `PDS_oracle`).

Semisynthetic continuous sweep (so the x-axis is continuous, not four points): take real
cells and DIAL the window by (i) subsampling depth m and (ii) scaling the true effect
magnitude (multiply `delta_v - delta_bar_g` by a factor), tracing `D_self/D_null` continuously
from ~1 down to JAK1-like values. Show the recovery metric is a single function of the window
and the four real genes land on it.

Interpretation and expected central result: when `D_self/D_null` ~ 1 (TP53/KRAS), even large
true alpha-differences yield indistinguishable PDS (tau ~ 0, high inversion) -> the benchmark
cannot recover model order; as the window widens, tau -> 1, inversions -> 0. The crossing
point defines a "benchmarkable" threshold.

Inputs: build/truth-split pseudobulks, `harness.PDS_cos`. Outputs:
`TBD(run:controlled_predictors.py -> results/controlled_recovery.csv, results/resolution_curve.csv)`.

### 1D. Analytical scaling law g(rho) and collapse (theory layer)

Purpose: explain WHY effect size is the predictive feature (the honest Fig 5a observation) and
put 1A/1B/1C on one axis; the inverse of g is the design calculator.

Derivation (state assumptions explicitly in Methods; derive the SCALING, calibrate the SHAPE):
- Pseudobulk sampling model: half-mean x_bar ~ N(mu_v, Sigma_v / m). Then
  `E||x_bar_a - x_bar_b||^2 = ||mu_a - mu_b||^2 + tr(Sigma_a)/m + tr(Sigma_b)/m`.
- With per-profile noise `eta^2 = tr(Sigma)/m`: `D_self^2 = 2 eta^2`,
  `D_null^2 = ||delta||^2 + 2 eta^2`, so `(D_self/D_null)^2 = 1/(1 + S)`,
  `S = ||delta||^2 / (2 eta^2)`. Sanity check with canonical ratios: JAK1 r=0.21 -> S ~ 21.7,
  TP53 r=0.96 -> S ~ 0.085, GATA1 r=0.88 -> S ~ 0.29 (self-consistent).
- PDS discriminates v from other variants, so the relevant signal is the between-variant
  spread `Delta_bar_v = RMS_u ||delta_v - delta_u||`. The score depends on the data only
  through the dimensionless group `rho = Delta_bar * sqrt(n) / (2 s)` and the competitor count
  n_test: `E[PDS] = g(rho; n_test)`, g monotone, g(0) = 0.5.

Honest boundary: we DERIVE the scaling (rho ~ Delta * sqrt(n) / s, dependence only through rho
and n_test) and CALIBRATE g by the 1C semisynthetic sweep, then confirm the four real genes and
the gene-level atlases fall on g. We do NOT claim a closed form for tie-aware cosine PDS. This
also explains why raw depth n alone is a weak predictor: Delta varies ~100x across genes while
sqrt(n) does not, so effect size dominates (matching the honest predictor finding that adding
log-cells barely helps).

Output: `TBD(run:floor_law_fit.py -> results/floor_law.csv)`.

---

## 6. Phase 2: defensive robustness (reviewer-survival)

### 2A. Metric-family generalization
Recompute the floor quantity (`D_self/D_null`, `PDS_oracle`, or rankable fraction) under
several distributional metrics so a reviewer cannot dismiss it as energy-distance small-sample
behavior: energy distance, MMD (RBF), Wasserstein/Sinkhorn, pseudobulk cosine/L2, and a
classifier two-sample test (cross-validated AUROC of a logistic/RF classifier trained to
separate variant-A cells from variant-B cells at the single-cell level). The classifier test
is the most model-free "is the signal there at all" probe and pairs naturally with 1B. Message:
the floor is metric-family-invariant; each family has its own noise floor.
Output: `TBD(run:metric_family_floor.py -> results/metric_family_floor.csv)`.

### 2B. Pilot-to-full prospective validation
Per perturbation: use only 25-50 pilot cells -> estimate pilot effect size and noise ->
predict rankability (and, with 1D, required depth) at 100/200/400 cells -> validate on held-out
independent cells. Report a calibration curve and predicted-vs-observed depth threshold rather
than only AUROC. This directly answers the "LODO is retrospective / same-source" critique and
is what upgrades "triage" toward a genuine (still-bounded) power tool.
Output: `TBD(run:pilot_to_full.py -> results/pilot_to_full.csv)`.

### 2C. Replicate-aware noise (conditional on data)
Random split-half shares batch, library and editing, so it is an optimistic lower bound on true
measurement noise. IF a dataset carries replicate/batch structure (verify first for Ursu
TP53/KRAS, GATA1 base-editing, JAK1 scSNV-seq), compare `D_replicate` vs `D_split_half` and show
split-half is optimistic, the replicate-aware floor is stricter, and the ordering is preserved.
This can only strengthen the null. If no dataset has usable replicates, do not fake it; instead
state in text that split-half is an optimistic lower bound.
Output: `TBD(run:replicate_noise.py -> results/replicate_noise.csv)` or a text-only caveat.

### 2D. Continuous calibrated reliability score
Replace the hard `S > W` headline with a continuous, calibratable quantity, e.g.
`Reliability_v = P(D_null > D_self + margin)` estimated by bootstrap, with user-selectable
operating points (exploratory < 0.5, uncertain 0.5-0.9, benchmarkable > 0.9). Keep the binary
label as a chosen operating point, not the definition. Do this alongside 1C (it is cosmetic
without the resolution-curve validation).

---

## 7. Phase 3: reach (do not gate the paper on these)

- 3A. One independent variant dataset that breaks the gene x technology x cell-type confound
  (a weak-effect high-depth or a replicate-carrying dataset would be ideal). RUNX1 is BLOCKED
  (needs raw SRA reprocessing or author-provided per-cell genotypes). Strongest form if it
  lands: freeze/pre-register the rankability rule, then validate on the new dataset as external
  validation without re-running all models.
- 3B. Deeper gene-level atlas panel: 7-10 Perturb-seq datasets (CRISPRi/a/KO, varying depth,
  replicate structure) from the scPerturb portal; test whether measurement-reliability ->
  model-ranking-stability holds across datasets, generalizing the framework beyond allele level.
- 3C. Title/positioning upgrade. DONE 2026-07-30, ahead of 1C/1D: the title is now
  "Measurement and model limits of allele-specific single-cell perturbation prediction",
  and the abstract, introduction, Results, Discussion and figure captions were repositioned
  so that allele-specific prediction, not AllelePerturb, is the paper's subject. The new
  title states scope rather than a result, so it does not depend on 1C/1D landing; both
  limbs (measurement and model) are required, because attributing the failure to
  measurement alone contradicts the JAK1 model-limited regime.

---

## 8. Target main-figure specifications (text; author renders)

One-line message per figure, then the panels as text with axes and the traced result file.
Fig 4 is new (the keystone); Fig 3 gains the oracle and pairwise panels; Fig 5 becomes an
evidence-bearing workflow rather than a schematic.

| Fig | Message (one sentence) | Panels (axes / content) | Traces to |
|-----|------------------------|--------------------------|-----------|
| 1 | Allele-level prediction is a distinct, harder task than gene-level | a task schematic gene vs allele; b dataset table (gene, aa, #variants, tech, cell type); c per-variant depth (median cells/variant, log y); d representation and split schematic; e evaluation axes (direction vs discrimination) | dataset tables; `allele_perturb_bench.csv` |
| 2 | Direction recovery does not imply allele discrimination | a paired Pearson-delta (0.60-0.68) vs PDS (0.49-0.52) per predictor with 95% CI; b PDS-vs-Pearson-delta scatter, all predictors in one quadrant; c DE gradient (direction 78% > LFC 46% > overlap 28%); d per-gene dissociation (TP53 0.30, KRAS 0.18, GATA1 0.15, JAK1 0.02); e explicit direction-only baseline showing high correlation, chance PDS | `results_v4_10metrics.csv`; add direction-only run |
| 3 | Ground truth has a finite measurement window that bounds what is knowable | a split-half concept; b D_self/D_null per gene (0.96, 1.00, 0.88, 0.21 with CI); c NEW oracle PDS ceiling per gene and vs depth (1A); d NEW pairwise allele-allele resolvable fraction vs detectable-vs-WT fraction (1B); e depth curves (detection rate vs cells/half) | canonical ratios; `oracle_ceiling.csv`; `pairwise_resolvability.csv`; `split_half_power_curve.csv` |
| 4 (NEW) | The measurement window predicts whether a benchmark can recover the true model ordering | a controlled-predictor PDS vs alpha per gene/regime (flat for TP53/KRAS, monotone for JAK1); b benchmark-resolution curve: ranking-recovery (Kendall tau or 1-inversion) vs window (D_self/D_null or rho), real genes + semisynthetic sweep on one curve, "benchmarkable" threshold marked; c analytical g(rho) collapse across genes + atlases (1D); d metric-family robustness (floor across energy/MMD/Wasserstein/classifier-2-sample); e gene-level atlas extension of reliability -> stability | `controlled_recovery.csv`; `resolution_curve.csv`; `floor_law.csv`; `metric_family_floor.csv` |
| 5 | A power-aware, prospective workflow triages what is worth benchmarking | a honest effect-size AUROC per dataset (Replogle 0.97, Norman 0.93, Adamson 0.93, GATA1 0.95, JAK1 1.00; TP53/KRAS not evaluable); b NEW pilot-to-full calibration (predicted vs observed rankability / depth, 2B); c depth as an effect-size-conditioned design variable; d power-aware triage workflow with the reliability operating points (2D) | `rankability_predictor_honest.csv`; `pilot_to_full.csv` |

Supplementary / Extended Data: full 18x3 head grid; alternative rankability criteria; replicate
noise (2C) if available; per-metric floor tables; the derivation of g(rho) with assumptions.

---

## 9. New tables (draft skeletons; numbers filled from runs, never guessed)

### Table A. Per-gene measurement-window summary (main-text candidate)
Consolidates detection, identification, ceiling and model performance in one honest view.

| Gene | D_self/D_null (95% CI) | PDS_oracle (ceiling, 95% CI) | % variants identifiable | % allele-pairs resolvable | % rankable (native) | best-model PDS |
|------|------------------------|------------------------------|-------------------------|---------------------------|---------------------|----------------|
| TP53 | 0.96 (0.95-0.98) | 0.485 (0.445-0.525) | 0.0 | 0.0 | 0.0 | ~0.49 |
| KRAS | 1.00 (0.99-1.02) | 0.500 (0.440-0.566) | 0.0 | 0.0 | 0.0 | ~0.49 |
| GATA1| 0.88 (0.86-0.90) | 0.572 (0.549-0.595) | 0.0 | 3.8 | 2.4 | ~0.49 |
| JAK1 | 0.21 (0.12-0.33) | 0.792 (0.742-0.838) | 15.4 | 78.5 | 90 | ~0.42-0.51 |

Executed 2026-07-10 (oracle_ceiling.csv, pairwise_resolvability.csv). Reading: TP53/KRAS
oracle ceilings sit AT chance (CI crosses 0.5) -> pure measurement limitation, no model can
help. JAK1 ceiling is 0.79 while its best model is <=0.51 -> measurement is adequate there and
the residual gap is genuinely computational (a real, honest split of the two limitations).

### Table B. Controlled-predictor benchmark resolution (Fig 4 support)
Bootstrap over held-out variants; recovery is a sharp function of the ceiling (window).
Executed 2026-07-10 (controlled_recovery.csv).

| Regime (ceiling = oracle window) | P(correct order) | P(winner correct) |
|----------------------------------|------------------|-------------------|
| ceiling <= 0.52 (TP53 all, KRAS most) | 0.27 | 0.38 |
| ceiling 0.52-0.56 (GATA1 modest window) | 0.95 | 0.97 |
| ceiling 0.65-0.80 | 1.00 | 1.00 |
| ceiling > 0.80 (JAK1) | 0.97 | 0.97 |

Second axis: recovery also falls when the number of held-out variants n_var is small
(GATA1 m=250, n=73 -> 0.80; JAK1 m=150, n=6 -> 0.92), matching g(rho; n_test).

### Table C. Metric-family floor (SI or Fig 4d)
| Metric family | TP53 | KRAS | GATA1 | JAK1 | floor persists? |
|---------------|------|------|-------|------|-----------------|
| energy distance | 0.96 | 1.00 | 0.88 | 0.21 | yes (baseline) |
| MMD (RBF) | TBD(2A) | TBD(2A) | TBD(2A) | TBD(2A) | TBD |
| Wasserstein/Sinkhorn | TBD(2A) | TBD(2A) | TBD(2A) | TBD(2A) | TBD |
| classifier 2-sample AUROC | TBD(2A) | TBD(2A) | TBD(2A) | TBD(2A) | TBD |

---

## 10. Risk register

| Risk | Mitigation |
|------|------------|
| Semisynthetic collapse looks tautological (noise and rho same-sourced) | Noise from real held-out cells; estimate s, Delta from a split independent of scoring; use real delta directions; validate g on the four real genes without fitting to them |
| Controlled-predictor alpha=1 trivially self-matches -> inflated ceiling | Build delta from a separate build-half; alpha=1 ceiling equals PDS_oracle by construction, not 1.0 |
| CLT/Gaussian breaks at very low n (the regime that matters) | 1C sweep covers low n empirically; g is calibrated not extrapolated; report where the closed-form scaling deviates |
| cosine vs Euclidean geometry mismatch in the derivation | Derive scaling in Euclidean, then show the same collapse empirically for tie-aware cosine PDS; flag as a stated subtlety |
| Framework ambition dilutes the clean null | All Phase 1 experiments are metrology bounds; they strengthen or bound the null and never produce a model-performance claim |
| Over-claiming a power calculator | Keep "triage and diagnostic" until 2B calibration exists; the g(rho) inverse is presented as bounded guidance, not exact cell counts |
| gene x tech x cell-type confound read as biology | Phase 0 caveat #4; frame JAK1 as a favorable regime, not stronger biology |
| Canonical-scorer drift across runs | P1 locks one scorer before any Phase 1 run |

---

## 11. Sequencing and dependencies

1. P1, P2 (scorer + S/W definition) -> unblocks all of Phase 1.
2. Phase 0 text (independent, ship anytime).
3. 1A -> 1C (ceiling defines the honest alpha=1 anchor) -> 1D (rho axis) ; 1B in parallel with 1A.
4. Fig 3 update (b + new 1A, 1B panels) and Fig 4 assembly (1C, 1D) after Phase 1 runs.
5. Phase 2 (2A, 2B, 2D) after Fig 4 exists; 2C only if replicate data verified.
6. Phase 3 opportunistic; 3C title only after 1C/1D are in the manuscript.

Definition of done for "structurally NM": Fig 4 exists and shows a monotone benchmark-resolution
curve (ranking recovery rises with the measurement window), the four real genes fall on the
analytical g(rho), and the abstract/Discussion claims are scoped to what these experiments
actually establish.
