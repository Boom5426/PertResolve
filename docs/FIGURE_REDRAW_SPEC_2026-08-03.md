# Figure redraw specification: printed Figures 5 and 6

> **Status 2026-08-03.** Printed Figures 3 and 4 were resized to these targets
> and are done: both pass the panel gate and are in the manuscript. Figures 5
> and 6 were resized too, but ten of their fourteen panels then failed on
> overlapping labels, which is design work rather than arithmetic. The
> remaining work is listed under "What is still open" at the end.

Written 2026-08-03 for the panel designer.

## The problem in one paragraph

The `Figure3_refined_v2`, `Figure4_refined_v2`, `Figure5_panels_refined` and
`Figure6_panels_refined` deliveries each state that they "deliberately do not contain a full-figure
composite". The composite was therefore never designed, and the panel canvases
were chosen without a total-height budget. Placed at their own declared sizes
and packed optimally, the four figures come to 277.4, 241.8, 379.0 and 503.6 mm
tall against budgets of 187.6, 174.2, 160.9 and 183.0 mm. Scaling them down at placement is
not available: every panel already bottoms out at 5.0 to 5.6 pt type, so any
reduction at placement time drops it below the 5 pt Nature Methods floor. The
panels have to be redrawn on smaller canvases.

Nothing else is wrong with them. All 30 panels render cleanly, embed one font
family, and read their numbers from the committed result tables; those fixes are
already applied in the repository and are described at the end so they are not
undone by the redraw.

## Budgets, and where they come from

A figure that is legal at 183 mm wide can still be unplaceable, because the
float carries its caption as well. The captions differ in length, so each figure
has its own budget. These were measured, not estimated, by substituting a blank
183 mm by H page for the real figure and bisecting on whether LaTeX reports
"Float too large":

| file | printed as | max composite height | currently | must lose |
|---|---|---|---|---|
| `fig1.pdf` | Figure 1 | 165.4 mm | 163.0 | done |
| `fig2.pdf` | Figure 2 | 187.6 mm | 159.0 | done |
| `fig4.pdf` | **Figure 3** | **187.6 mm** | 277.4 | **89.8 mm** |
| `fig3.pdf` | **Figure 4** | **174.2 mm** | 241.8 | **67.6 mm** |
| `fig5.pdf` | **Figure 5** | **160.9 mm** | 379.0 | **218.1 mm** |
| `fig6.pdf` | **Figure 6** | **183.0 mm** | 503.6 | **320.6 mm** |

Reproduce with `python scripts/figures/measure_figure_budget.py`. Re-run it if a
caption changes length: a longer caption shrinks its figure's budget.

Note the file-to-number mapping. `fig3/` and `fig4/` are swapped relative to the
printed numbers, as `manuscript/figures/README.md` records. The deliveries were
named by printed number, so they were installed crosswise.

## Hard constraints

1. **Width 183 mm.** A row is panels plus a 4.5 mm gutter between them, so the
   panels in one row must sum to at most 178.5 mm for two panels, 174.0 for
   three.
2. **Type floor 5.0 pt, absolute.** Type is set in points on a canvas declared
   in millimetres, so shrinking the canvas does not shrink the type. Text will
   collide instead. Where a target is much smaller than the current canvas,
   expect to remove content rather than compress it.
3. **The canvas is exact.** Write the declared millimetre canvas, not a tight
   bounding box: `nm_style.save(fig, path, exact=True, formats=("pdf", "png",
   "svg", "tiff"))`. The compositor places every panel at the width it was drawn
   at and applies no scale factor, so a trimmed canvas silently changes the
   on-page type size.
4. **Rows are 2.6 mm apart**, which is what the 8 pt panel letter occupies.
   The band above row 1 costs another 2.6 mm, so a figure's height is
   `sum(row heights) + 2.6 x number_of_rows`. An earlier version of this
   document omitted that leading band and every target came out 2.6 mm
   optimistic; the numbers below include it.
5. **Panel letters must read a, b, c ... in order** across rows, so a row can
   only contain panels that are already adjacent in the sequence.

## Targets

Widths that are unchanged from the current delivery are marked `=`. The row
heights are the binding numbers; the widths only matter where a row is currently
too wide to fit at all.

### printed Figure 3, repository `fig4/`, budget 187.6 mm  (DONE)

Achieved 187.4 mm. g and h ended at 51 mm rather than 52 to pay for the leading
letter band.

| row | panels | target size (mm) | row height | now |
|---|---|---|---|---|
| 1 | a `fig4a_pds_splits`, b `fig4b_direction_splits` | 76 = x **44**, 76 = x **44** | 44 | 54 |
| 2 | c `fig4c_empirical_null`, d `fig4d_distance_invariance` | 89 = x **37**, 62 = x **37** | 37 | 45 |
| 3 | e `fig4e_feature_spaces`, f `fig4f_gene_split_dissociation` | 76 = x **45**, 76 = x **45** | 45 | 55 |
| 4 | g `fig4g_model_harness`, h `fig4h_interface_audit` | **94** x **52**, **84** x **52** | 52 | 64 |

Total 44 + 37 + 45 + 52 + 3 x 2.6 = **185.8 mm** of 187.6.

Row 4 is the one that also has to narrow: g and h are currently 100 and 89 mm,
which is 189 mm before any gutter and cannot share a row at any height.

### printed Figure 4, repository `fig3/`, budget 174.2 mm  (DONE)

Achieved 172.4 mm at row heights 40 / 43 / 43 / 36, one millimetre under each
target below for the same reason.

| row | panels | target size (mm) | row height | now |
|---|---|---|---|---|
| 1 | a `fig3a_window_definition`, b `fig3b_representative_windows` | 55 = x **41**, 68 = x **41** | 41 | 58 |
| 2 | c `fig3c_window_ratio`, d `fig3d_sibling_identification` | 70 = x **44**, 89 = x **44** | 44 | 62 |
| 3 | e `fig3e_classifier_control`, f `fig3f_effect_depth_landscape` | 89 = x **44**, 70 = x **44** | 44 | 62 |
| 4 | g `fig3g_depth_rankability`, h `fig3h_native_rankability` | 76 = x **37**, 72 = x **37** | 37 | 52 |

Total 41 + 44 + 44 + 37 + 3 x 2.6 = **173.8 mm** of 174.2.

Every width already fits; this figure needs height only. d and e are currently
36 mm tall and are being asked to grow to 44, which is free room for their
labels; the cost falls on a, b, c and f.

### printed Figure 5, repository `fig5/`, budget 160.9 mm

This is the severe one. Five of its seven panels are currently 98 to 112 mm
wide, so no two of them can share a row and each takes a full row of its own.
Getting to four rows means bringing six panels under 90 mm as well as cutting
their heights roughly in half.

| row | panels | target size (mm) | row height | now |
|---|---|---|---|---|
| 1 | a `fig5a_replicate_ceiling`, b `fig5b_ceiling_vs_models` | **90** x **35**, **88** x **35** | 35 | 52, 64 |
| 2 | c `fig5c_phase_map`, d `fig5d_predictor_family` | 70 = x **35**, **88** x **35** | 35 | 64, 50 |
| 3 | e `fig5e_rank_recovery_concept`, f `fig5f_external_benchmarks` | **90** x **37**, **88** x **37** | 37 | 50, 68 |
| 4 | g `fig5g_recovery_support` | 112 = x **45** | 45 | 82 |

Total 35 + 35 + 37 + 45 + 3 x 2.6 = **159.8 mm** of 160.9.

If the content cannot survive that compression, say so rather than shrinking the
type: the fallback is to keep a, b and c in the main figure (121.2 mm at current
sizes, which fits today) and move d to g to Extended Data. That is a scientific
call for the author, not a layout one, so raise it rather than deciding it.


### printed Figure 6, repository `fig6/`, budget 183.0 mm

The worst of the four, and the only one where **narrowing is unavoidable**. All
seven panels are 100 to 140 mm wide, so the cheapest possible pair is
100 + 4.5 + 100 = 204.5 mm and no two of them can share a row at any height.
Seven rows of one comes to 503.6 mm. Dropping panels does not rescue it either:
even a, b and c alone, each on its own row, is 203.2 mm, already over budget.
Six of the seven have to come under 90 mm so they can pair.

| row | panels | target size (mm) | row height | now |
|---|---|---|---|---|
| 1 | a `fig6a_pilot_pipeline`, b `fig6b_retrospective_roc` | **89** x **41**, **89** x **41** | 41 | 58, 68 |
| 2 | c `fig6c_dataset_auroc`, d `fig6d_prospective_validation` | **89** x **49**, **89** x **49** | 49 | 72, 82 |
| 3 | e `fig6e_effect_depth_landscape`, f `fig6f_paired_depthmatch` | **89** x **44**, **89** x **44** | 44 | 74, 66 |
| 4 | g `fig6g_probability_workflow` | 140 = x **41** | 41 | 68 |

Total 41 + 49 + 44 + 41 + 3 x 2.6 = **182.8 mm** of 183.0.

d is the hardest single case: 132 x 82 mm down to 89 x 49, a 60% area cut. If
the prospective-validation panel cannot carry its calibration curve and its AUC
summary at that size, splitting it into two 89 mm panels on separate rows is
better than shrinking the type, and the row budget above has no slack for that,
so raise it rather than deciding it.

## Defects fixed in the repository: please do not reintroduce

The v2 deliveries shared a set of defects that are now fixed in
`manuscript/figures/`. Start from the repository copies, not from the delivery.

1. **Font stack.** Each `figN_common.py` restated the whole style with
   `["Arial", "Helvetica", "DejaVu Sans"]` and `mathtext.rm = "Arial"`. On a
   machine without Arial that renders every panel in DejaVu Sans, silently and
   without error. Worse, the Arial that Debian/Ubuntu's
   `ttf-mscorefonts-installer` provides has no U+0302, so `\hat{...}` takes its
   accent from STIXGeneral and puts a second typeface in the panel. Style now
   comes from `nm_style.resolve_sans()`, which picks a family by glyph coverage
   and raises when none qualifies. `fig4a_pds_splits.py` additionally overrode
   the stack inside the panel; that is removed too.
2. **Two symbols have no glyph in any Arial-metric sans**, and mathtext
   substitutes STIXGeneral for them without warning: `\cup` (U+222A) and `\ll`
   (U+226A). They are spelled out instead. Check any new symbol against
   `nm_style.REQUIRED_CODEPOINTS`.
3. **Private data copies.** Each delivery shipped its own `data/` directory. A
   panel that reads its own copy keeps plotting the old numbers after the
   committed table is regenerated. Reads now resolve by name through
   `nm_style.data_path()`, which searches the repository and raises rather than
   falling back. All 19 tables the deliveries used resolve to committed files.
4. **`pearson_delta_bootstrap_summary.csv` had no generator.** The recipe was
   recovered (`numpy.random.default_rng(0)`, 2,000 resamples, percentile 2.5 and
   97.5, one shared generator across methods in a fixed order) and is now
   `scripts/figures/make_pearson_bootstrap_summary.py`, which reproduces the
   shipped values exactly and verifies them with `--check`. Any new derived
   table needs the same treatment.
5. **Two panels overran their canvas** by 0.22 and 0.38 mm
   (`fig4c_empirical_null`, `fig4g_model_harness`); both are fixed.
6. **The delivered `benchmark_resolution.py` and `allele_resolution.py` were
   older than the repository's**, still carrying hard-coded `/data/boom/...`
   paths that the repository version had already replaced with `--atlas-dir`
   and `--out`. They were not copied back.
7. **`fig4b_selfnull_dist.csv`** (200 seeds by 2 variants) exists only in the
   delivery; it is now committed at `results/canonical/fig4b_selfnull_dist.csv`.
   Its gene-level ratios agree with the committed `bootstrap_CIs.json`
   (TP53 P222P 0.927 against a gene mean of 0.965; JAK1 R108Q 0.164 against
   0.210, CI 0.116 to 0.327), but **its generator is still only on the compute
   server** and should be synced for full traceability.
8. **Figure 6 is the one delivery that got data provenance right**, and it is
   the model to copy: a `data/raw` tree that is byte-identical to the committed
   tables, a `data/derived` tree, and a `prepare_fig6_data.py` that regenerates
   the second from the first. It was verified by deleting the derived tables and
   re-running the preparation in the compute environment: five of nine come back
   byte-identical and the other four agree to 2e-4, which is logistic-regression
   solver noise across library versions (the ROC geometry itself, `fpr` and
   `tpr`, is bit-identical; only the fitted `threshold` moves). The tables are
   committed at `results/fig6_derived/` and the preparation script at
   `scripts/figures/prepare_fig6_data.py`.

## What is still open

Figures 5 and 6 were resized to the targets above and every panel renders, but
ten of the fourteen then failed the gate on **overlapping labels**: at 35 to
49 mm tall these panels carry more prose than fits, so an annotation lands on the
plot or on its neighbour. This is the part that needs redrawing rather than
rescaling. The gate reports each one; as of 2026-08-03 they are:

| panel | collision |
|---|---|
| `fig5a_replicate_ceiling` | title and footer both run past the 90 mm width |
| `fig5b_ceiling_vs_models` | axis label over "The maximum is a selected point" |
| `fig5c_phase_map` | the two regime annotations overlap completely |
| `fig5e_rank_recovery_concept` | title runs past the 90 mm width |
| `fig5f_external_benchmarks` | axis label over two separate footnotes |
| `fig5g_recovery_support` | "depth m range" over the gene and depth labels |
| `fig6a_pilot_pipeline` | "perturbation + reference" runs off the left edge |
| `fig6b_retrospective_roc` | "JAK1" and "18 / 2" over the single-class footnote |
| `fig6d_prospective_validation` | "Calibration" over "50-cell training-free SNR" |
| `fig6e_effect_depth_landscape` | gene counts over the subtitle |

`fig5a`, `fig5e` and `fig5g` are the hardest: they are schematics drawn in
absolute coordinates, so the artwork rescales with the canvas while the type does
not, and the prose has to be shortened or re-placed rather than moved. The rest
are annotation placement.

## Verification

```bash
python manuscript/figures/check_panels.py 5 6
```

It fails on: text outside the declared canvas or within 0.4 mm of a side edge,
two labels overlapping by more than half of the smaller one, a PDF that is not
its declared millimetre canvas, any type below 5 pt, a row wider than 183 mm, a
figure taller than its budget, and any panel embedding more than one font family.
The side margin is not zero because matplotlib's text bounding box excludes an
italic's slant overhang: a label reported at 93.79 mm on a 94.0 mm canvas was
still being cut. Update the
declared sizes in its `SPECS` table when the canvases change, since that table
is what the composite arithmetic is checked against.

Then build the composite with

```bash
python manuscript/figures/compose.py fig4 fig4_assemble a,b c,d e,f g,h
```

which writes the `.tex`, runs `lualatex`, and prints whether the result fits.
