# Figure 5: assembly guide

**One-line message:** the measurement window sets an achievable ceiling and governs
benchmark validity: a benchmark cannot rank models reliably when the ground-truth
discrimination ceiling is itself near chance. Source of truth:
`manuscript/latex/AllelePerturb_manuscript.tex`. House palette via shared
`manuscript/figures/nm_style.py`.

## Panel status (current pass, all panels built and placed)

Every panel is a matplotlib script that emits a vector PDF at (or within a few percent
of) its placement width in `fig5_assemble.tex`, so on-page type stays 5.5-7 pt and the
composite embeds zero raster. The three schematics are no longer AI raster imports; the
`*_prompt.md` files are kept only as the content spec they were drawn from.

| Panel | Content | Placed | File / source |
|-------|---------|--------|---------------|
| a | oracle-ceiling definition (2nd measurement) | 64.00 mm | `fig5a_oracle_def.py` (schematic, 64 x 41.1 mm; no numbers) |
| b | per-gene oracle ceiling vs best model (dumbbell) | 56.82 mm | `fig5b_oracle_vs_model.py` : `results/canonical/oracle_ceiling.csv` + `unified_results5.csv` |
| c | measurement-limited vs computation-gap phase map | 44.19 mm | `fig5c_phasemap.py` : same two canonical tables as b |
| d | synthetic-predictor construction (α interpolation) | 85.25 mm | `fig5d_synthetic.py` (schematic, 85 x 36 mm; no numbers) |
| e | low- vs high-resolution leaderboard example | 85.25 mm | `fig5e_leaderboard.py` (schematic, 85 x 36 mm; no numbers) |
| f | benchmark resolution vs oracle ceiling, 9 datasets | 85.25 mm | `fig5f_resolution.py` : `results/benchmark_resolution/summary.csv` |
| g | **keystone**: benchmark-validity transition | 85.25 mm | `fig5g_validity.py` : `results/canonical/controlled_recovery.csv` |

Legacy raster panels `fig5a_oracle_def.pdf`, `fig5d_synthetic.pdf`, `fig5e_leaderboard.pdf` are superseded
and no longer referenced by the assembly.

## Layout and encoding contract (2026-07-27 Nature Methods compliance pass)

- Page is exactly **183 x 150 mm** (`border=0mm`; previously 185 x 152 mm because the
  1 mm standalone border was added outside the 183 mm bounding box).
- Margins 3.5 mm, **every gutter 5.5 mm**. Row 1 (a|b|c) shares one top edge at
  y = 146.3 and one baseline (all three panels are 41.04-41.10 mm tall by construction).
  Row 2 (d|e) and row 3 (f|g) are exact-height pairs.
- Row 3 (f, g) is the largest row (85.25 x 57.9 mm each): the ceiling result (b, c) and
  the validity transition (f, g) are the message, the schematics support them.
- **Panel letters**: lowercase bold at `\fontsize{8}{9}`, i.e. 8 pt on the page. The old
  `\bfseries\large` rendered ~12 pt.
- **One legend rule**: no gene legend anywhere. Gene identity is direct-labelled in b
  (coloured y-tick labels), c, f and g. Panel g previously carried two legends (a gene
  key plus the metric key); the gene key is gone, the metric key stays.
- **Gene hues are reserved**: TP53 #5185C0, KRAS #E99D4E, GATA1 #8281B9, JAK1 #55966B
  mean those genes and nothing else. Schematics a, d, e therefore use only the house
  slate ramp (`FEATURE_COLORS` #9AA7B3 -> #5F6B76 -> #2E3742), GREY and INK. The muted
  red that panel e used for the scrambled ordering is removed (there is no red in the
  house palette); the scramble now reads from crossing dashed leaders plus a broken
  lightness order, which also survives greyscale.
- **Marker shape encodes one variable only.** Panel f used shape for dataset class while
  the neighbouring panel g uses shape for the metric; f is now all circles and dataset
  class is carried by direct labels. Fill state follows one rule in b and g: filled =
  the primary/limiting quantity (oracle ceiling; full-order recovery), open = the
  secondary one (best model; winner selection).
- **f and g are a matched pair**: identical x-label wording, identical x-limits
  (0.46-0.98), y-limits, ticks and regime bands. The bands are labelled once, in f.
- d and e are a matched pair: identical canvas, identical caption strip geometry, the
  same 5.5 pt grey sub-headers / 5.8 pt labels / 5.6 pt caption type ladder, and e
  re-uses d's P1..P5 ramp so the two panels thread together.
- Every panel PDF is generated at (or within 0.3% of) its placement width, so the point
  sizes written in the scripts are the on-page point sizes: median 6.1-6.2 pt, minimum
  4.7 pt (two mathtext subscripts in d), nothing above 7.5 pt, zero embedded raster.

## 5f verified numbers (self-contained, committed)

`benchmark_resolution/summary.csv`, (oracle_ceiling, resolution_P_recover_order),
reproduces SI Table 5 exactly. Values as committed in that CSV (re-read 2026-07-27):
TP53 0.492→0.014, KRAS 0.493→0.018, GATA1 0.520→0.512, sci-Plex 0.602→0.772,
VCC 0.731→1.00, Replogle 0.783→1.00, Adamson 0.836→1.00, Norman 0.848→1.00,
JAK1 0.887→0.99. (An earlier revision of this README listed KRAS as 0.499→0.036; the
committed CSV, and therefore the panel, say 0.493→0.018. Nothing was changed in the
data, only this note corrected.) Below a ceiling near 0.5 (TP53/KRAS) resolution
collapses; above ~0.6 (GATA1 partial, then JAK1 + gene-/drug-level screens) it reaches ~1. The x-axis here is
the ≥100-cell / 50-per-half benchmark-resolution oracle (SI Table 5 depth), NOT the
native-depth per-gene oracle of the main-text window table (that is panel 5b).

## Two different ceilings live in this figure

- Panels **b** and **c** use the **native-depth per-gene oracle ceiling** with 95% CIs
  (`results/canonical/oracle_ceiling.csv`: TP53 0.485, KRAS 0.500, GATA1 0.572,
  JAK1 0.792), against the best in-house head from `unified_results5.csv`.
- Panel **f** uses the **≥100-cell / 50-per-half benchmark-resolution ceiling**
  (SI Table 5 depth), so JAK1 sits at 0.887 there and 0.792 in b. Both axis labels say
  so; do not merge the two.
- Panel **g** uses the per-depth controlled sweep (`controlled_recovery.csv`), whose
  GATA1 ceilings span 0.517-0.549. The caption's "clears ~0.57 (GATA1)" refers to the
  native-depth ceiling of panel b, not to g's x-axis; flagged for the text pass.

## Notes

- All panels import the shared `nm_style`; vector PDF (0 raster) + 600 dpi PNG.
- Schematics a/d/e carry no results by design (no ceiling value, no model or dataset
  name); the quantitative statements live in b, c, f, g.
- Re-render the composite with `pdflatex fig5_assemble.tex`, then copy
  `fig5_assemble.pdf` over `fig5_composite.pdf` and `../../latex/figures/fig5.pdf`.
- Do NOT put split-half ratio, external model audit, feature comparison, rankability
  predictor or the pilot workflow in Fig 5 (those are Fig 3/4/6).

## Open items for the text pass (figure was NOT changed to hide these)

1. **b and c plot the same eight numbers** (per-gene oracle ceiling and per-gene best
   in-house PDS), in two forms. c earns its place only through the 2-D regime taxonomy
   and the y = x reference; b now no longer repeats the regime words, so the taxonomy
   lives in exactly one panel. If a panel has to be cut, c is the merge candidate: its
   content is b plus the diagonal. Not done here because the caption fixes the a-g
   roster and the caption is out of scope for this pass.
2. **Caption g says "clears ~0.57 (GATA1)"** but panel g's x-axis is the per-depth
   controlled sweep, whose GATA1 ceilings are 0.517-0.549; 0.572 is the native-depth
   ceiling of panel b. The number is correct for b, wrong for g. Text-side fix.
3. **Caption b says the best model's "interval also overlapp[s] 0.50"**, but panel b
   draws a CI only for the oracle, not for the model. Either add the model CI (needs a
   committed per-gene CI for the best head) or reword. Not invented here.
4. Panel c keeps an empty upper-left area. That emptiness is the result (no model
   reaches its ceiling). The y-axis is now tightened to 0.46-0.72 while x still spans
   every ceiling (0.46-0.82), so equal aspect was dropped; the y = x diagonal is drawn
   explicitly and the regime shading is defined from it, so above/below the diagonal
   stays unambiguous. y cannot be tightened to the plotted models (max 0.54) without
   deleting the high-ceiling stretch of the diagonal that the caption's
   "benchmark-solvable regime" refers to.
5. **Panels f and g are different quantities and must stay labelled as such.** f is one
   point per DATASET at a matched shallow depth (nine benchmarks, Methods "External
   benchmark calibration": 50 cells per half, perturbations with >=100 cells); g is one
   point per GENE x SEQUENCING DEPTH for the four allele genes (15 cells, depth 25-250)
   and carries two metrics. They disagree by construction where they look comparable
   (GATA1: resolution 0.512 in f versus P_correct_order 0.989 at depth 50 in g), so the
   panel titles and x-axis labels now state the unit of analysis. Neither panel is a
   duplicate; do not merge them.
