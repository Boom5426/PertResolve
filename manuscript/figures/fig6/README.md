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
| a | pilot → rankability concept | **matplotlib schematic, drawn at final size** | `fig6a_pipeline.py` → `fig6a_pipeline.pdf` (60 × 45 mm, composite scale 1.0). Supersedes the AI raster-model draft `fig6a_pipeline.pdf` (423 mm canvas placed at 60 mm, i.e. ~2.4 pt on-page type); `fig6a_prompt.md` is kept as the content brief only |
| b | LODO ROC curves | **data-direct done** | `fig6b_roc.py` — reconstructed from `second_probe_rankability_table.csv`; AUCs match `rankability_predictor_honest.csv` exactly |
| c | dataset-level AUROC forest | **data-direct done** | `fig6c_auroc.py` — `rankability_predictor_honest.csv` (effect_size) |
| d | prospective pilot-to-full validation (quant. half) | **data-direct done** | `fig6d_prospective.py` — `pilot_validation/pilot_validation_summary.json` (T100) + `pilot_validation/README.md` (SNR, null) |
| e | effect-size-conditioned design landscape | **data-direct done** | `fig6e_design.py` — `native_rankability` (effect_size × n_cells × rankable) |
| f | native vs depth-matched (n=50) external atlases | **data-direct done** | `fig6f_depthmatch.py` — `canonical_numbers.json` |
| g | power-aware reporting protocol, **decision half only** | **matplotlib schematic, drawn at final size** | `fig6g_workflow.py` -> `fig6g_workflow.pdf` (168 x 52 mm, placed at 176 mm). Supersedes the external vector art `fig6g_workflow.pdf`; `fig6g_prompt.md` is kept as the content brief only |

Run any panel: `python fig6<x>_*.py` (imports shared `nm_style` + `fig6_data`); vector PDF
(0 raster) + 600 dpi PNG.

## a / g scope split (2026-07-27 rework)

Panels a and g used to tell the same story: g's first three boxes ("Pilot screen",
"Estimate measurement properties", "Predict rankability") were exactly panel a's three
steps, and g took roughly the bottom half of the page to restate them in external vector
art whose font, rounded boxes, red and green arrows and warning triangle clashed with the
matplotlib panels. The two are now sequential and share one visual language:

* **a = what you measure and predict.** It ends at the predictor OUTPUT (a predicted
  probability of clearing the split-half floor) and carries the Results boundary
  condition verbatim: detection-level triage, perturbation versus reference, necessary
  but not sufficient for allele identification. It no longer states any action.
* **g = what you then do.** It starts at the decision node, with an entry chip that
  *names* panel a instead of redrawing it, and owns the Yes / No branches, the re-pilot
  feedback and the stratified report.

Both are drawn in matplotlib at their exact placement size, in the same box idiom and
type sizes. Palette: neutral slate / grey ramp plus house HOTSPOT amber (#E69F00) for the
single decision node. No red (there is none in the house palette) and no gene hue: green
#55966B stays JAK1's colour in b, c and e, so Yes / No are separated by position, label
and box weight and the panel survives greyscale.

## Figure-wide encoding rules

* One colour language, house palette only (`fig6_data.DATASET_COLORS`). A hue means a
  *gene*: TP53 #5185C0, KRAS #E99D4E, GATA1 #8281B9, JAK1 #55966B, in b, c and e alike.
  The external gene-level atlases are not genes, so they take the house neutral slate ramp
  (`nm_style.FEATURE_COLORS`) instead of hues of their own: Replogle #2E3742, Norman
  #5F6B76, Adamson #9AA7B3, dark to light in the order in which they rank across c/d/f.
* One dataset key only: a bare colour-to-name list in **b** (no title, no numbers), plus
  directly coloured axis labels in **c** and **f**. The AUROC values and their intervals
  are stated once, in c; the TP53 / KRAS single-class note is carried once, in c.
* Greyscale redundancy in **b**: the five ROC curves also carry one dash pattern each
  (`fig6b_roc.DASHES`), none of them the chance diagonal's, and the key's handle is long
  enough to show a full dash period, so the panel does not rely on hue alone.
* One mark, one fill rule in **d/e/f**: circle everywhere; filled = the primary condition
  (learned predictor / rankable / native depth), open = the comparison condition
  (training-free predictor / un-rankable / depth-matched). The square marker formerly in d
  was the only competing shape encoding and has been removed.
* Reference-line idiom shared by **c** and **d**: dashed grey = null or chance, dotted ink
  = the mean.
* Direct labels wherever the categories are spatially stable (d's two predictors, e's four
  genes, f's two depth conditions), so only b and e still carry a legend, one each.

**Legibility rule for schematics:** draw them in matplotlib at their *final* placement
size (composite scale ≈ 1.0), as `fig6a_pipeline.py` does, so in-panel type stays 5-7 pt
on the page. `fig6a_pipeline.py` therefore turns off the shared `savefig.bbox="tight"`
so the exported PDF is exactly the 60 × 45 mm box used in `fig6_assemble.tex`.

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
  cutoff or fabricated contour (positioned as a triage landscape). The y axis is the
  number of cells **scored at the deepest split-half rung** (n_work = 150, so n_cells =
  300), not the raw sequencing depth; 341 of the 464 perturbations sit exactly on that
  cap, which the axis label and an explicit hairline cap guide now state. Readability was
  fixed by marker weight (small, semi-transparent open markers for un-rankable, larger
  opaque markers for rankable) and by moving the gene names into a header band above the
  cap; no point was jittered and no value was changed.
- **6f depth-match** (`canonical_numbers.json`): un-rankable % native→matched-50 =
  Replogle 55.3→54.9, Adamson 14.6→38.5, Norman 3.4→11.2. Depth normalization raises
  Adamson/Norman toward Replogle, but Replogle stays high (dataset-specific structure).

## Framing guards (kept in the panels)

Position as a **triage tool, not a universal cell-number calculator**: wording is
"predicted rankability / relative triage / effect-size-conditioned design / expected
measurement regime", never "required cell number / guaranteed rankability / universally
sufficient depth". No split-half ratio, external model audit, feature comparison, or
oracle content in Fig 6 (those are Fig 3/4/5).

Panel a carries the boundary condition stated in the Results: "Detection-level triage:
perturbation vs reference, necessary but not sufficient for allele identification". Panel
g stays at "measurement floor" / "measurable perturbations" wording and never claims
allele-level identification, so the two panels are consistent with the text.
