# Figure 1, assembly guide

**One-line message:** AllelePerturb defines perturbation prediction at protein-coding
variant resolution (nested solvability, coverage, sampling regime, features, evaluation,
splits); a definition figure, **no results**.

Source of truth is the LaTeX manuscript `manuscript/latex/AllelePerturb_manuscript.tex`.
The older `manuscript/figures/` composites (`fig1_final.*`, `fig1_benchmark_overview.*`)
are **superseded**; do not reuse them.

## Composite

`fig1_composite.pdf` (assembled by `fig1_assemble.tex`, build: `lualatex fig1_assemble.tex`).
Layout: row 1 **a | b**, row 2 **c | d (schematic + PCA)**, row 3 **e | f**, with bold
lowercase 8 pt panel letters **a-f**; page 185 x 165.0 mm (bounding box 183 x 163.0 mm).
Every placed file is vector; the composite contains zero embedded raster
(`pdfimages -list fig1_composite.pdf | tail -n +3 | wc -l` -> 0).

Each row is flush to the 183 mm page on both edges and sums to exactly 183 mm:

| Row | Panels | Widths (mm) | Gutter |
|---|---|---|---|
| 1 | a, b | 113.5 + 66.0 | 3.5 |
| 2 | c, d schematic, d PCA | 58.0 + 69.0 + 47.0 | 4.5 |
| 3 | e, f | 78.0 + 100.0 | 5.0 |

Every panel is placed at exactly the width it was drawn at, so its internal 5-7 pt type is
5-7 pt on the page; no placement carries a corrective scale factor. Panel scripts write
their declared canvas through `nm_style.save(exact=True)`.

**Do not change a panel's width without re-running its script.** The panels are drawn in
absolute millimetres with type in points, so scaling one at placement time shrinks its type
below the 5 pt Nature Methods floor. Panel a is 113.5 mm rather than the 122.5 mm it was
first drawn at for exactly this reason: at 122.5 mm row 1 came to 188.5 mm before any
gutter, and neither a nor b could absorb the difference by scaling because both bottom out
at 5.0 pt.

## Quality gate

`python check_fig1_panels.py` is the gate. It re-runs every panel in memory and fails on:

- text outside the declared canvas (this is how the clipped legend footer in b and the
  overrunning "own observed" label in e were found);
- a figure whose size is not its declared millimetre canvas;
- any text below 5 pt;
- a composite row wider than 183 mm.

It also writes `../nm_font_resolved.tex`, which records the family the panels actually
resolved to so `nm_fonts.tex` stamps the panel letters in the same one. Run it after
editing any panel, then rebuild the composite.

## Panel inventory (each panel is its own file)

| Panel | Content | Type | Script | Output | Data source |
|-------|---------|------|--------|--------|-------------|
| a | nested detection / identification / prediction, with the measurement-limited vs model-limited verdict map | matplotlib schematic | `fig1a_nested_solvability.py` | `fig1a_nested_solvability.*` | none |
| b | variant positions on 4 proteins, tracks and key in one file | **data-direct** | `fig1b_variant_tracks.py` | `fig1b_coverage.*` | `data/allele_perturb_bench_v2.csv` |
| c | per-variant cell depth | **data-direct** | `fig1c_depth.py` | `fig1c_depth.*` | `allele_perturb_bench_v2.csv` (`n_cells`) |
| d | variant feature concept + shared feature space | matplotlib schematic + **data panel** | `fig1d_theta_schematic.py`, `fig1d_theta_pca.py` | `fig1d_theta_schematic.*`, `fig1d_theta_pca.*` | theta columns of `allele_perturb_bench_v2.csv` |
| e | evaluation decomposition, PDS primary | matplotlib schematic | `fig1e_eval_axes.py` | `fig1e_eval_axes.*` | none |
| f | six generalization splits | **data-direct** | `fig1f_splits.py` | `fig1f_splits.*` | `split*_role` columns of `allele_perturb_bench_v2.csv` |

Each panel writes `.pdf` (deliverable), `.png` (600 dpi preview), `.svg` (editable) and
`.tiff` (LZW, 600 dpi).

**The data source is `allele_perturb_bench_v2.csv`, not `allele_perturb_bench.csv`.**
The two differ in `is_hotspot` on 182 of the 470 rows. In the default table that column is
a leaked continuous score that is not even an indicator (it sums to -6.16 over TP53 and
+11.27 over KRAS); v2 replaces it with the external-annotation definition and keeps the old
values as `is_hotspot_leaked_OLD`. Methods states the external definition
("defined exclusively from external prior knowledge, with no reference to the single-cell
outcome"), and the reported evaluation grid `results/results_v4_exttheta.csv` is computed
on the de-leaked theta, so drawing Fig. 1 from the default table put the figure at odds
with both. Pass `corrected=True` to `nm_style.load_bench`.

## Style (Nature Methods)

One sans family per figure, resolved by `nm_style.resolve_sans()`, 5-7 pt, 0.5 pt axes,
no top/right spines, vector PDF with editable text. Single-column 89 mm / double 183 mm.
Resolution is by glyph coverage, not just by whether a family is installed: the Arial from
Debian/Ubuntu `ttf-mscorefonts-installer` has no U+0302, so every `\hat{...}` in this
figure would take its accent from STIXGeneral and put a third typeface in the panel. On
this machine the stack therefore resolves to Liberation Sans, which is Arial-metric and
covers the whole set.

Palette is the house set in `manuscript/figures/nm_style.py` (`GENE_COLORS`), shared by
every figure in the manuscript, **not** Okabe-Ito:
TP53 `#5185C0` (blue), KRAS `#E99D4E` (orange), GATA1 `#8281B9` (purple),
JAK1 `#55966B` (green); hotspot/held-out accent `#E69F00`, grey `#7A7A7A`,
light grey `#D9D9D9`, ink `#1A1A1A`. One gene keeps one colour across all panels.
(An earlier version of this README documented an Okabe-Ito set that the panels never
used; corrected 2026-07-27.)

### Legibility rule (enforced)

Every panel is drawn at its exact placement width, so its placement scale is 1.00 and its
source type reaches the page unchanged. `check_fig1_panels.py` fails on any text below
5 pt and on any panel whose PDF is not its declared millimetre canvas, which is what makes
the scale-1.00 claim checkable rather than asserted.

## Numbers that appear (all trace to the bench CSV)

Real protein-coding variants **470** (TP53 98, KRAS 92, GATA1 254, JAK1 26); median
cells/variant **929 / 1000 / 354 / 104**; protein lengths 393 / 189 / 413 / 1154 aa;
technologies Perturb-seq (TP53, KRAS) · base editing (GATA1) · scSNV-seq (JAK1).

> **RESOLVED (2026-07-27):** panel c prints the GATA1 median as **354** (WT-excluded,
> 470-variant convention, reproduced by `python fig1c_depth.py`) and the manuscript Results
> now also says "354 for GATA1". Figure and text agree.

JAK1 domain labels in panel b are the short module names **JH2** / **JH1** (pseudokinase /
kinase), matching the naming used in this script's docstring and Supplementary Table 1; the
long forms "JH2 pseudokinase" / "JH1 kinase" overlapped each other at 5.5 pt.

> **RESOLVED, variant count convention (2026-07-11): 470.** The committed CSV has 472
> rows but 2 are WT reference rows (KRAS WT 644 cells; GATA1 WT 38,276 cells), so there
> are **470** real coding variants (TP53 98, KRAS 92, GATA1 254, JAK1 26). These panels
> use the 470-variant (WT-excluded) convention via `nm_style.load_bench(exclude_wt=True)`,
> which also removes a spurious 38 k-cell "variant" from the GATA1 depth panel. The
> manuscript text and README were updated to 470 (KRAS 92, GATA1 254), and total-cell
> statements now note the 321,043 includes wild-type/control cells. Set
> `exclude_wt=False` only to reproduce the old 472/93/255 numbers.

## Figure legend

The caption lives in `manuscript/latex/AllelePerturb_manuscript.tex` only. A second copy
here would drift from it; the panels changed on 2026-08-03 and the copy that used to sit
in this file still described the superseded panel a.

## Legend economy (Nature Methods)

Figure 1 carries exactly **two** keys, one per encoded variable, and **no gene legend**:

* **b**, inside `fig1b_coverage.pdf`: marker shape = variant consequence (missense /
  nonsense / synonymous), amber fill = external hotspot annotation. The "colour = gene"
  half was removed because every track is already direct-labelled with its gene name in
  the gene colour.
* **f**: light grey = training variants, steel blue (`nm_style.HELDOUT`) = held-out. Amber
  is **not** used here; it means external hotspot annotation everywhere in this figure and
  must not also mean "held out".

The gene to colour mapping is carried by **direct labels only**: the gene-coloured track
titles in b and the gene-coloured tick labels in c. The PCA in d has no legend of its own;
it direct-labels each gene cluster and additionally varies the marker (circle, square,
triangle, diamond), because `GENE_COLORS` is not separable in greyscale or under
deuteranopia (see the `nm_style` docstring), so no panel in this figure encodes gene by
hue alone.
