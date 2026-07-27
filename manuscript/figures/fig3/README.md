# Figure 3 — assembly guide

**One-line message:** split-half analysis reveals an allele-resolution measurement
floor (within-variant noise vs variant signal); TP53/KRAS/GATA1 sit near the floor,
JAK1 retains a wide window. Source of truth: `manuscript/latex/AllelePerturb_manuscript.tex`.

Palette: house set (TP53 `#5185C0`, KRAS `#E99D4E`, GATA1 `#8281B9`, JAK1 `#55966B`)
via the shared `manuscript/figures/nm_style.py`.

## Panel status (Nature Methods compliance pass)

Composite letters are the caption's; file stems keep their historical suffix.

| Letter | File | Unique question | Source |
|--------|------|-----------------|--------|
| a | `fig3a_window_def.py` | What *is* the detection window? | schematic, no data |
| b | `fig3b_windows.py` | What does a closed vs an open window look like for one variant? | `results/canonical/fig3b_selfnull_dist.csv` (200-seed split-half; TP53 P222P R=0.93, JAK1 R108Q R=0.16) |
| c | `fig3c_ratio.py` | Where does every variant of every gene sit relative to the floor? (hero) | `second_probe_rankability_table.csv` + `bootstrap_CIs.json` |
| d | `fig3g_detect_ident.py` | Can sibling alleles be told apart by distance? | `results/canonical/pairwise_resolvability.csv` |
| e | `fig3h_classifier.py` | Can a supervised single-cell classifier tell them apart? | `results/canonical/classifier_two_sample.csv` |
| f | `fig3d_landscape.py` | Is it depth or effect size that sets the window? | native-depth rows (effect_size, D_null, n_cells) |
| g | `fig3e_titration.py` | Does adding cells close the window? | `split_half_power_curve.csv` (PCA-50) |
| h | `fig3f_rankable.py` | What fraction of variants clears the detection bar at native depth? | `unrankable_canonical.json` |

Run any panel: `python fig3<x>_*.py` (all import the shared `nm_style` and, for data
panels, `fig3_data.py` or `remote_data.py`). Each writes a vector PDF (0 embedded
raster) + 600 dpi PNG preview.

## Layout contract (do not break when editing a panel)

- **Every panel is drawn natively at its placed width**, so the composite scale is
  exactly 1.00 and on-page type stays in the 5-7 pt house band. If you change a
  panel's `W_MM`/`H_MM`, re-measure its PDF page size and update the matching
  `\includegraphics[width=...]` in `fig3_assemble.tex`.
- Rows are **bottom-aligned** (common baseline); left/right margins 3 mm; gutters
  6.2 mm (row 1), 8.0 mm (row 2), 7.6 mm (row 3). Composite 183 x 167 mm.
- Panel letters: lowercase bold, `\fontsize{8}{9}`, `anchor=north west`.
- **One key per panel, at most.** The gene colour key is never repeated: genes are
  direct-labelled (coloured bold tick labels in c/d/e/h, cloud labels in f,
  end-of-curve labels in g). The only keys in the figure are the two style keys
  (d: pale vs full-strength bar tint; e: open vs filled marker) and f's marker-area key.
- **Open/filled means exactly one thing in the figure, and only in e**: open =
  detection (vs wild type), filled = identification (sibling vs sibling), as the
  caption states. Panel d's two series are both sibling-level identification
  measures, so d codes them by tint of the gene hue (pale = the looser pair-level
  measure, full = the stricter variant-level measure) and its key shows filled
  swatches only. Do not reintroduce an outline bar in d: it sits 8 mm from e's key
  and would make open/filled carry two meanings in one row.
- **D_self / D_null colour coding is shared by a and b**: D_null is grey in both;
  D_self is the emphasised quantity (ink in the monochrome schematic a, gene hue in
  the data panel b). Do not invert this in a.
- Zero-valued bars (d, h) sit on a light 0-100 % track so an honest 0 % reads as a
  measured zero on the full range rather than as empty canvas.
- **d and h must not look like the same object**: d is paired horizontal bars,
  h is vertical columns with bootstrap CIs. They are both four-gene 0-100 %
  summaries and were being read as one panel twice.

## Numbers (all verified to trace to committed files)

- D_self/D_null (mean, 95% CI): TP53 0.96 (0.95–0.98), KRAS 1.00 (0.99–1.02),
  GATA1 0.88 (0.86–0.90), JAK1 0.21 (0.12–0.33). Native-depth aggregation reproduces
  `canonical_numbers.json` exactly.
- Rankable at native depth: TP53 0%, KRAS 0%, GATA1 2.4%, JAK1 90%
  (n = 98 / 92 / 254 / 20 evaluated perturbations).
- Depth titration (PCA-50, fraction detectable): JAK1 = 1.00 at all rungs; TP53
  0.51→0.69 (no 300 rung), KRAS 0.57→0.74, GATA1 ~0.51–0.57 (plateau < 75%).
- Effect size (median): TP53 1.41, KRAS 1.11, GATA1 2.79, JAK1 13.04.

## Notes / decisions

- **External gene-level atlases moved out of Fig 3.** Per the agreed plan, Replogle
  (55.3%), Norman (3.4%), Adamson (14.6%) un-rankable belong in Fig 6 / Extended Data,
  not here; Fig 3 stays allele-focused. `unrankable_canonical.json` still holds those
  values for the ED/Fig 6 panel. The manuscript currently shows them in Fig 3d — that
  text/structure edit is pending sign-off.
- The formerly "held" panels b/d/e now load from committed `results/canonical/`
  tables via `remote_data.py`; nothing in Fig 3 depends on an uncommitted file.
- 3a is programmatic (not AI) to stay vector and palette-consistent. It is drawn at
  its composite placement size (48 x 55 mm, scale 1.00), so its type lands at 6 pt
  on the page; the retired AI `fig3a_window_def.pdf` was 635 mm wide and rendered at 1.5 pt.
  It is deliberately monochrome (ink / grey): the four gene hues are reserved for
  the data panels, so a neutral schematic cannot collide with a gene colour.
- **d/h repeat-reading (addressed without a caption edit):** panels d and h carry
  near-identical numbers (d resolvable 0/0/3.8/78.5 %, h rankable 0/0/2.4/90 %) and
  used to be the same object twice (four-gene horizontal bars, same order, same
  0-100 % axis, same track). They measure different things (identification vs
  detection) and the Results cite them separately, so the presentation, not the
  data, was the problem: h is now drawn as vertical columns and each panel names its
  own quantity on its value axis ("Sibling-allele identification (%)" in d,
  "Variants rankable at native depth (%)" in h). Merging them into one three-series
  panel would still be cleaner but needs a caption relettering, so it stays for
  sign-off. `fig3f_rankable.py` `H_MM` is now 65.85 (not 58.5) purely so that, after
  the orientation flip, its placed height keeps row 3 bottom-aligned.
- **Retracted "manuscript number to check" note (2026-07-27):** this file previously
  claimed the Results text's GATA1 D_self/D_null CI (0.85-0.90) disagreed with
  `bootstrap_CIs.json`. It does not. `dself_dnull_ci["GATA1"]` is
  lo = 0.85477, hi = 0.90295, which rounds to 0.85-0.90 exactly as the prose states.
  No manuscript edit is needed; do not "fix" the lower bound to 0.86.
