# Figure 2 — assembly guide

**One-line message:** variant-level prediction reveals a direction-ranking dissociation:
every predictor recovers the perturbation direction but none ranks the correct allele.
Source of truth: `manuscript/latex/AllelePerturb_manuscript.tex`. House palette + shared
`nm_style.py`; data via `remote_data.py`.

## Panel status

| Panel | Content | Type | File / source |
|-------|---------|------|---------------|
| a | unified-evaluation schematic | **redrawn, vector** | `fig2a_evaluation.py` (matplotlib, drawn at its final 59 x 63 mm placement size; no numbers) |
| b | PDS forest (20 in-house predictors) | **data done** | `fig2b_pds_forest.py` - `results/canonical/definitive_summary.csv`; owns the shared label column + shared key |
| c | Pearson-δ forest | **data done** | `fig2c_pearson_forest.py` - `results/results_v4_exttheta.csv`; reuses `fig2b.row_order()` and `fig2b.interval()` |
| d | PDS-vs-Pearson scatter | **data done** | `fig2d_scatter.py` — PDS from definitive_summary, Pearson from exttheta |
| e | DE fidelity gradient | **data done** | `fig2e_de_gradient.py` — exttheta; drawn with reproducible 65/0.25/0.15, manuscript text updated |
| f | per-gene dissociation gap | **data done** | `fig2f_gene_gap.py` — exttheta |
| g | representative TP53 variants | **data done** | `fig2g_pervariant.py` — exttheta (Ridge-esm TP53) |

Feature-space colour (2b/2c/2d) uses the shared `S.FEATURE_COLORS` slate ramp
(theta = light slate, ESM = mid slate, ESM+theta = dark slate, reference = grey), one
consistent mapping across the three panels and deliberately distinct from the gene hues,
which are now used only where the variable really is the gene (2f, 2g).

## Nature Methods compliance pass (composite v3)

Page is now **183 x 170 mm** (`border=0mm`; the previous build shipped 185 mm wide).
Every panel PDF is drawn at its final placement size, so every `\includegraphics`
width equals the panel's own page size and the placement scale is 1.00.

1. **One encoding per variable.** 2b and 2c previously separated the four feature-space
   categories with marker shapes (b) and dash patterns (c). Both now use one circle with
   a white edge plus the slate ramp, and `fig2c` imports `interval()` and the `LW/MS/CAP`
   constants from `fig2b`, so an interval is drawn identically in the two panels. 2d's
   bespoke diamond/square reference markers were replaced by the same grey circle.
2. **One key, one label column.** 2b draws the only feature-space key in the figure; 2c
   and 2d draw none (2d's colours are direct-labelled by its four callouts). 2c reuses
   `fig2b.row_order()` and sets no y tick labels, so the 20 method names appear once.
   2b and 2c share the identical axes band (bottom 8 mm, height 47.5 mm), so their rows
   line up across the gutter.
3. **Palette discipline.** 2a no longer borrows TP53 blue and KRAS orange for its two
   evaluation branches (those hues mean genes in 2f/2g and in Figs 1/3/6); the schematic
   is ink + grey with one `S.HOTSPOT` highlight on the "own variant" candidate. 2e's blue
   ramp (#8FB8DE/#5185C0/#2C5A8F, a TP53-hue reuse) became a single neutral grey, because
   colour was not encoding anything there.
4. **Axes and dead space.** 2d dropped the decorative quadrant shading and the unexplained
   y = 0.30 divider, and its limits were tightened from 0.42-0.58 to 0.482-0.523. 2e's
   x label is now self-explanatory ("Agreement with measured response (1 = perfect)").
   2g was rebuilt as two real stacked axes with real y-axis labels on a common 0-1 scale
   (the old top band stretched 0.70-0.85 across a strip, exaggerating the spread of the
   quantity the panel calls stable) and re-proportioned to 120.5 x 44 mm so it shares the
   bottom row with 2f instead of floating alone.
5. **Type.** Panel letters are lowercase bold at an explicit 8 pt (`\fontsize{8}{9}`),
   replacing `\bfseries\large` (~12 pt). Mathtext is pinned to Liberation Sans, so the
   composite embeds no DejaVu Sans and every glyph is Arial-metric. On-page type is
   5.5-7 pt; `pdfimages -list` reports 0 embedded rasters.

No plotted value changed in this pass.

## Legibility rework (composite v2)

The original `fig2a_evaluation.pdf` was a 1152 pt AI raster-style schematic placed at 48 mm
(scale 0.12), so its labels landed at ~2 pt on the page. It is superseded by
`fig2a_evaluation.py`, a matplotlib schematic drawn at the final placement size
(scale ~1.0, type 5.4-6.0 pt). `fig2a_evaluation.pdf` and `fig2a_prompt.md` are kept only as
the provenance of the panel's content.

Also in this pass: 2b/2c tick labels shortened to head + feature space
("KNN ESM+theta") and set at 6 pt in full-contrast ink (the reference rows keep
italics instead of grey), 2c's legend moved to the reserved strip above the axes
and its annotation off the label column, 2d's callouts relabelled to the same
short form, and 2g redrawn at 99 mm native width so it is placed at scale 1.0
instead of 1.42. Composite page is now 183 x 174 mm. No plotted value changed.

## Verified numbers

- **2b** PDS (multi-seed, `definitive_summary.csv`): 18 heads + Gene-mean + WT-null, range
  0.487-0.517, every 95% CI crosses 0.50 (external SOTA excluded -> Fig 4g).
- **2c** Pearson-δ (`results_v4_exttheta.csv`): per-method mean 0.555-0.647, all CIs > 0.
  Matches the manuscript's current Results wording ("Pearson-$\delta$ ranged from 0.55 to
  0.65 across predictors"). Panel c shows the 19 non-null predictors; WT-null is Fig 2b/2d
  only, so the caption phrase "the same predictors" covers 19, not 20, rows.
- **2d**: all predictors cluster in the "direction without ranking" region (PDS ~0.5, Pearson ~0.6); WT-null at (0.50, 0).
- **2f** gap = mean(Pearson) - mean(PDS): TP53 0.31, KRAS 0.19, GATA1 0.13, JAK1 -0.04 (~ manuscript 0.30/0.18/0.15/0.02; JAK1 ~ 0).
- **2g** Ridge-esm TP53: Pearson stable 0.735-0.821, PDS 0.000-0.948 across the 10 shown
  held-out variants. Note the Results sentence quotes PDS "0.00 to 0.98", which is the range
  over **all 84** Ridge-esm TP53 rows (max 0.9794), not over the 10 representatives plotted
  (max 0.9485); the Pearson range it quotes (0.74 to 0.82) is the plotted subset. Not changed
  here, since the panel is described as representative; flagged for the author.

## RESOLVED: panel 2e (DE gradient)

The manuscript DE gradient 78/46/28 reproduced from **no committed or remote file**
(results_v4 gives ~65/0.25/0.15; no method/aggregation on the server reaches 78/46/28).
Per the author's decision (option A), 2e is drawn with the reproducible in-house exttheta
values **direction 65% / DE-LFC Spearman 0.25 / DE overlap 15%**, and the manuscript text
was corrected in two places (Results + Fig 2 caption). Same decreasing gradient, lower
absolute values.
