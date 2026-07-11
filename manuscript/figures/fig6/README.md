# Figure 6 — assembly guide

**One-line message:** rankability prediction enables power-aware evaluation — pilot the
assay, estimate measurement resolution, predict rankability, size or redesign the
experiment, and benchmark models only where the ground truth clears the measurement
floor. Source of truth: `manuscript/latex/AllelePerturb_manuscript.tex`. House palette +
shared `manuscript/figures/nm_style.py`; shared loaders in `fig6_data.py`.

**Figure 6 is the most locally complete figure: 5/7 panels are data-direct and done.**

## Panel status

| Panel | Content | Type | File / source |
|-------|---------|------|---------------|
| a | pilot → rankability concept | AI schematic | `fig6a_prompt.md` |
| b | LODO ROC curves | **data-direct done** | `fig6b_roc.py` — reconstructed from `second_probe_rankability_table.csv`; AUCs match `rankability_predictor_honest.csv` exactly |
| c | dataset-level AUROC forest | **data-direct done** | `fig6c_auroc.py` — `rankability_predictor_honest.csv` (effect_size) |
| d | prospective pilot-to-full validation (quant. half) | **data-direct done** | `fig6d_prospective.py` — `pilot_validation/pilot_validation_summary.json` (T100) + `pilot_validation/README.md` (SNR, null) |
| e | effect-size-conditioned design landscape | **data-direct done** | `fig6e_design.py` — `native_rankability` (effect_size × n_cells × rankable) |
| f | native vs depth-matched (n=50) external atlases | **data-direct done** | `fig6f_depthmatch.py` — `canonical_numbers.json` |
| g | **hero**: power-aware workflow (decision flow) | AI schematic | `fig6g_prompt.md` |

Run any panel: `python fig6<x>_*.py` (imports shared `nm_style` + `fig6_data`); vector PDF
(0 raster) + 600 dpi PNG. The two AI panels (6a concept, 6g workflow): paste prompt into a
raster model, overlay the exact text in Illustrator (6g is text-heavy — overlay, don't
trust generated text).

## Verified numbers (all trace to committed files)

- **6b/6c LODO AUROC** (effect_size, `rankability_predictor_honest.csv`): Replogle 0.97,
  Norman 0.93, Adamson 0.93, GATA1 0.95, JAK1 1.00; mean **0.96**. TP53/KRAS
  `evaluable=False` (all un-rankable, single-class) — omitted from the plot, noted only.
  6b reconstructs the LODO logistic (log₁₀ effect size) and its AUCs match the committed
  values to 0.00.
- **6d prospective** (`pilot_validation_summary.json` T100): learned effect-size
  Adamson 0.85, Norman 0.94, Replogle 0.99 (mean **0.93**); training-free pilot-SNR
  0.84/0.93/0.98 and label-permutation null **0.49** from the pilot README; near-diagonal
  calibration from the JSON `calibration` bins. **Red line:** the pooled
  `mechanistic_auroc` 0.978 is base-rate-inflated (README says do-not-report) and is
  deliberately excluded; only the 3 balanced datasets are shown.
- **6e design landscape** (`native_rankability`, 464 allele-gene perturbations): rankable
  fraction TP53 0/98, KRAS 0/92, GATA1 6/254, JAK1 18/20; high-effect variants rankable
  even at ~60 cells, low-effect variants un-rankable at 300. No universal cell-number
  cutoff or fabricated contour (positioned as a triage landscape).
- **6f depth-match** (`canonical_numbers.json`): un-rankable % native→matched-50 =
  Replogle 55.3→54.9, Adamson 14.6→38.5, Norman 3.4→11.2. Depth normalization raises
  Adamson/Norman toward Replogle, but Replogle stays high (dataset-specific structure).

## Framing guards (kept in the panels)

Position as a **triage tool, not a universal cell-number calculator**: wording is
"predicted rankability / relative triage / effect-size-conditioned design / expected
measurement regime", never "required cell number / guaranteed rankability / universally
sufficient depth". No split-half ratio, external model audit, feature comparison, or
oracle content in Fig 6 (those are Fig 3/4/5).
