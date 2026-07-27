# Figure 4: assembly guide

**One-line message:** allele-ranking failure persists across evaluation choices and
model interfaces (splits, metrics, representations, model class): an
alternative-explanation-elimination figure. Source of truth:
`manuscript/latex/AllelePerturb_manuscript.tex`. House palette via shared
`manuscript/figures/nm_style.py`.

## Nature Methods compliance pass (2026-07-27, second pass)

No plotted value changed. Structural/encoding changes only:

- **Panels a and b merged into one drawing, `fig4ab_splits.py` -> `fig4ab_splits.pdf`**
  (the old `fig4a_pds_splits.py` / `fig4b_pearson_splits.py` are removed). The two
  heatmaps show the same 18 heads x 5 splits, so the row-label column and the
  feature-space brackets are drawn ONCE and serve both maps. The letters `a` and `b`
  are stamped over the two blocks by `fig4_assemble.tex`.
- **One colormap logic for the a/b pair.** Previously a used a purple-green diverging
  map and b a blue sequential map for the same kind of quantity. Both now use ONE
  sequential lightness ramp (house neutrals -> dark slate) on ONE shared 0.20-0.80
  scale with ONE shared colourbar; the PDS chance level 0.50 is marked on the bar.
  Greyscale-safe (value maps monotonically to lightness).
- **Panel e**: the twin y axes were removed (both axes in fact spanned an identical
  0.30-0.75, so the second axis carried no information and made series-to-axis
  assignment ambiguous). Both series now share one y axis. The non-house terracotta
  red for "Direction" is gone; the Fig. 4 metric encoding is now direction =
  ink filled circle, PDS = grey open square, and the legend is replaced by direct
  labels.
- **Panel f**: the four-entry gene legend is replaced by direct gene labels next to
  each gene's points; the "y = x" tag was dropped (the in-panel note and the caption
  already state what crossing the diagonal means).
- **Panel h**: green/red status fills replaced by a monotone lightness ladder from
  house tokens (dark slate / HOTSPOT amber for the single caveat / LIGHT_GREY /
  near-white; luminance 105/162/219/240). Meaning is carried by the vector glyph
  SHAPE (tick / wave / cross / dot), so the matrix survives greyscale and CVD.
- **Palette discipline across the figure**: gene hues now appear ONLY in the one
  gene-resolved panel (f). Panels c, d, e, g were recoloured from the TP53 blue /
  terracotta accents to house neutrals (GREY / LIGHT_GREY / INK), so no hue does two
  jobs anywhere in Fig. 4.
- **Composite**: panel letters are now lowercase bold at an explicit 8 pt
  (`\fontsize{8}{9}\selectfont\bfseries`) instead of `\large` (~12 pt). Every panel
  is placed at scale 1.00; panels within a row share a common top edge and a common
  baseline to within 0.8 mm; gutters are equal within each row. Limits were tightened
  in c (0.35-0.60, no data clipped: observed range 0.357-0.589) and d (0.40-0.52,
  data range 0.441-0.504). Page is unchanged at 183 x 166 mm, 0 embedded raster.

## Legibility/label pass (2026-07-27, first pass)

All eight panels are drawn and composited (`fig4_assemble.tex` -> `fig4_composite.pdf`
-> `latex/figures/fig4.pdf`). The table below is the earlier data-sourcing status and is
kept for provenance; it no longer reflects "drawable now".

What changed in this pass (no plotted value was altered):

- **Metric name corrected.** `fig4d` y label and `fig4g` x label called PDS the
  "direction score"; PDS is the perturbation *discrimination* score. Both now match
  `fig2b_pds_forest.py` / `fig5b_oracle_vs_model.py`. Same correction in the `fig4d`
  and `fig4e` docstrings; `fig4f` x label was set to "PDS (cosine): allele discrimination" (shortened to
  "PDS (cosine)" in the second pass for cross-panel name consistency).
- **4a redrawn at final size** (native 60 x 49 mm, composite scale ~1.0). It now mirrors
  4b exactly: same row order (feature-space blocks x fixed head order), short head
  labels, block brackets, hand-built vector colourbar diverging about chance. The old
  version was an 18-row sorted heatmap with 4.4 pt in-cell numbers at scale 0.82.
- **4b, 4h** re-emitted at their placement size; 4h is now wide (109 x 52 mm) with
  horizontal column headers and aspect-corrected status glyphs.
- **4c** x label was clipped by its own tight bbox ("... per method x split x g"); it is
  now "Observed PDS (cosine) per combination".
- **4g** legend removed (redundant with the left-hand family labels, and it collided with
  the Gene-mean whisker); the "chance" tag and the panel message were moved clear of the
  whiskers.
- Every source font is now >= 5.0 pt and every panel sits at composite scale 0.97-1.09,
  so on-page type is 5.9-7.3 pt. Composite page is 183 x 166 mm (height reduced from 182).

## Panel status

| Panel | Content | Status | File / source |
|-------|---------|--------|---------------|
| c | empirical permutation-null calibration | **done** | `fig4c_null.py`, `permutation_null_pds.csv` (self-contained; verified null mean 0.500, 3.2% exceed) |
| h | interface-compatibility audit matrix | **done** | `fig4h_interface.py`, categorical, from `docs/REMOTE_INVENTORY.md` + Methods |
| a | PDS across 5 splits | **done** | `fig4ab_splits.py` (left block), `results/results_v4_exttheta.csv` |
| b | Pearson-δ across 5 splits | **done** | `fig4ab_splits.py` (right block), same source, shared rows and key |
| d | distance invariance (cosine/L1/L2) | **done** | `fig4d_distance.py`, `results/results_v4_exttheta.csv` (PDS_cos/L1/L2) |
| e | feature spaces (θ / ESM / ESM+θ) | **done** | `fig4e_features.py`, `results/results_v4_exttheta.csv` |
| f | gene × split dissociation | **done** | `fig4f_gene_split.py`, `results/results_v4_exttheta.csv` |
| g | **hero**: published external models under one harness | **done** | `fig4g_external_forest.py`, `results/canonical/definitive_summary.csv` via `remote_data.definitive()` |

**Scorer provenance (resolved).** The per-method forest (4g) reads the canonical
multi-seed/bootstrap summary in `results/canonical/definitive_summary.csv`, the same
source as Fig 2b, so the two figures quote one scorer. The per-split / per-gene /
per-feature breakdowns (4a, 4b, 4d, 4e, 4f) read the committed per-variant grid
`results/results_v4_exttheta.csv`, which is what reproduces the manuscript's breakdown
numbers (gap 0.31/0.19/0.13, GATA1 Low-depth PDS 0.19, 15 of 17). Breakdown PDS values
therefore sit at ~0.44-0.52 while the headline per-method values are 0.49-0.52; that is
a scorer/aggregation difference, not a discrepancy, and both are stated as such.

## 4c: verified numbers (self-contained, no remote dependency)

`permutation_null_pds.csv`, metric PDS_cos, n = 340 method×split×gene combinations:
permutation null mean **0.500** (per-combo means 0.495-0.504); **3.2%** exceed their
own 95th-percentile null (expected 5%); **33%** of combinations land exactly at 0.50
(tie-dominated small test sets). This closes the "maybe 0.50 is the wrong chance
baseline" objection.

## 4h: interface audit content (categorical facts)

Rows grouped as: **variant-conditionable** (scGen, scVIDR, Biolord, CellFlow ran to
completion and are scored; PerturbNet ran but its score is subspace-caveated) ·
**allele-blind by construction** (scGPT, GEARS, STATE condition on gene identity;
STATE was run per-gene but is still allele-blind) · **did not converge** (variant-CPA:
continuous θ input attempted, non-finite loss). Columns: continuous variant input /
unseen-allele conditioning / allele-specific output / ran to completion / allele-level
score defined. Gene-keyed resolution coverage = **4/470** variants (0.85%), updated
from 4/472 with the WT-row variant-count correction. Framing note: gene-keyed models
are an *interface incompatibility* (allele-level score undefined), NOT a performance
failure; keep that wording.

## Notes

- 4c and 4h import the shared `nm_style`; each writes a vector PDF (0 embedded raster)
  + 600 dpi PNG. 4h glyphs (check/cross/tilde) are drawn as shapes (Liberation Sans
  lacks ✓/✗), so they are font-independent.
- All eight panels regenerate from committed result files; none depends on a remote pull.
- Rebuild: run each `fig4*.py`, then
  `pdflatex -interaction=nonstopmode fig4_assemble.tex`, then copy `fig4_assemble.pdf`
  over `fig4_composite.pdf` and over `../../latex/figures/fig4.pdf`.
