# Figure 1 — assembly guide

**One-line message:** AllelePerturb defines perturbation prediction at protein-coding
variant resolution (task, coverage, sampling regime, features, evaluation) — a
definition figure, **no results**.

Source of truth is the LaTeX manuscript `manuscript/latex/AllelePerturb_manuscript.tex`.
The older `manuscript/figures/` composites (`fig1_final.*`, `fig1_benchmark_overview.*`)
are **superseded**; do not reuse them.

## Composite

`fig1_composite.pdf` (assembled by `fig1_assemble.tex`, build: `pdflatex fig1_assemble.tex`).
Layout: row 1 **a | b**, row 2 **c | d(+PCA inset) | e**, row 3 **f | g**, with bold panel
letters a-g; page 185 x 172 mm. Vector-preserving: the data panels (b/c/d-PCA) and Fig1f
stay vector; the AI raster panels are embedded at effective **460-791 dpi** at their placed
sizes (a 460, d 643, e 791, g 569 — all > 300, print-safe). To reposition anything, edit the
coordinates in `fig1_assemble.tex` and recompile. One easy Illustrator nudge: the d PCA inset
(`fig1d_theta_pca.pdf`) is dropped where the d schematic's "theta space" arrow points; move a
few mm if it clips the schematic label.

## Panel inventory (each panel is its own file)

| Panel | Content | Type | File(s) | Data source |
|-------|---------|------|---------|-------------|
| a | gene-level vs variant-level task | AI schematic | `fig1a_prompt.md` | — |
| b | variant positions on 4 proteins | **data-direct** | `fig1b_track_{TP53,KRAS,GATA1,JAK1}.pdf` + `fig1b_legend.pdf` | `data/allele_perturb_bench.csv`, `data/hotspot_external_definition.txt` |
| c | per-variant cell depth | **data-direct** | `fig1c_depth.pdf` | `data/allele_perturb_bench.csv` (`n_cells`) |
| d | variant feature concept | AI schematic + **data inset** | `fig1d_prompt.md` + `fig1d_theta_pca.pdf` | θ columns of the bench CSV |
| e | evaluation decomposition | AI schematic | `fig1e_prompt.md` | — |
| f | six generalization splits | AI schematic | `fig1f_prompt.md` | — |
| g | end-to-end workflow | AI schematic | `fig1g_prompt.md` | — |

Data-direct panels: `python fig1b_variant_tracks.py`, `python fig1c_depth.py`,
`python fig1d_theta_pca.py` (all import `nm_style.py`; PDF deliverable + 600 dpi PNG
preview). AI panels: paste the prompt into a raster image model, then overlay the exact
text in Illustrator (raster text is unreliable).

## Style (Nature Methods)

Helvetica/Arial (Liberation/Nimbus Sans stand-ins locally), 5–7 pt, 0.5 pt axes, no
top/right spines, vector PDF with editable text, colour-blind-safe Okabe-Ito genes:
TP53 `#0072B2`, KRAS `#D55E00`, GATA1 `#009E73`, JAK1 `#CC79A7`; hotspot `#E69F00`.
Single-column 89 mm / double 183 mm.

## Numbers that appear (all trace to the bench CSV)

Real protein-coding variants **470** (TP53 98, KRAS 92, GATA1 254, JAK1 26); median
cells/variant **929 / 1000 / 354 / 104**; protein lengths 393 / 189 / 413 / 1154 aa;
technologies Perturb-seq (TP53, KRAS) · base editing (GATA1) · scSNV-seq (JAK1).

> **RESOLVED — variant count convention (2026-07-11): 470.** The committed CSV has 472
> rows but 2 are WT reference rows (KRAS WT 644 cells; GATA1 WT 38,276 cells), so there
> are **470** real coding variants (TP53 98, KRAS 92, GATA1 254, JAK1 26). These panels
> use the 470-variant (WT-excluded) convention via `nm_style.load_bench(exclude_wt=True)`,
> which also removes a spurious 38 k-cell "variant" from the GATA1 depth panel. The
> manuscript text and README were updated to 470 (KRAS 92, GATA1 254), and total-cell
> statements now note the 321,043 includes wild-type/control cells. Set
> `exclude_wt=False` only to reproduce the old 472/93/255 numbers.

## Recommended figure legend (covers a–g)

**Figure 1 | AllelePerturb defines perturbation prediction at protein-coding variant
resolution. a**, Conceptual distinction between gene-level and allele-level perturbation
prediction. **b**, Benchmark coverage: variant positions along TP53, KRAS, GATA1 and
JAK1 with annotated domains; hotspot/pathogenic residues highlighted. **c**, Per-variant
cell depth spans distinct sampling regimes across genes (medians 929, 1000, 354, 104).
**d**, Each variant is represented by a 6-dimensional biophysical feature vector θ,
placing variants from all genes in a shared feature space. **e**, AllelePerturb-Eval
separates direction recovery, allele discrimination and differential-expression fidelity.
**f**, Six generalization splits. **g**, End-to-end benchmark workflow for held-out
protein-coding variant prediction.
