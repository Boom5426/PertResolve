# Figures

Figure-centric organization: each figure has its own directory holding the composite
(`figN_composite.pdf`), every individual panel as its own vector PDF plus a 600 dpi
PNG preview, the panel scripts, and a `README.md` recording the layout contract, the
data source of each panel and the decisions that must not be silently reverted.

## Directory to manuscript figure number

**Directory numbers are not manuscript figure numbers.** The Results order changed
after the files were named, and LaTeX numbers figures by order of appearance, so:

| Dir | File placed in the manuscript | Manuscript number | Subject |
|-----|-------------------------------|-------------------|---------|
| `fig1/` | `latex/figures/fig1.pdf` | **Figure 1** | task definition + benchmark coverage (no results) |
| `fig2/` | `latex/figures/fig2.pdf` | **Figure 2** | direction-vs-discrimination dissociation |
| `fig4/` | `latex/figures/fig4.pdf` | **Figure 3** | robustness across splits, metrics, representations, interfaces |
| `fig3/` | `latex/figures/fig3.pdf` | **Figure 4** | split-half measurement window, detection vs identification |
| `fig5/` | `latex/figures/fig5.pdf` | **Figure 5** | oracle ceiling and benchmark validity |
| `fig6/` | `latex/figures/fig6.pdf` | **Figure 6** | pilot rankability prediction and reporting protocol |
| `extended_data/` | `latex/figures/ED_fig1.pdf` | **Extended Data Fig. 1** | rankability verdict sensitivity to criterion choice |

`fig3/` and `fig4/` are therefore swapped relative to the printed numbers. Check the
`\includegraphics` line in the figure environment before editing a panel; do not infer
the target from the directory name.

## Building a composite

Each `figN/figN_assemble.tex` is a TikZ standalone that places the panel PDFs at their
final millimetre coordinates. Panels are drawn natively at their placed size, so every
placement scale is 1.00 and on-page type stays in the 5-7 pt house band.

```
cd figN && lualatex -interaction=nonstopmode figN_assemble.tex
cp figN_assemble.pdf figN_composite.pdf
cp figN_composite.pdf ../../latex/figures/figN.pdf
```

**lualatex, not pdflatex.** The panel letters must be set in the same sans family as
the panel type. `nm_fonts.tex`, shared by all six assemblies, resolves Arial ->
Helvetica -> Liberation Sans through `fontspec`, mirroring `nm_style.resolve_sans()`
on the matplotlib side. pdflatex cannot load a system OpenType/TrueType family and
will fail here by design; the previous `\usepackage{helvet}` route succeeded under
pdflatex but embedded NimbusSanL (a Helvetica clone) for the letters alongside
Liberation Sans (an Arial clone) for everything else, i.e. two typefaces per figure.
Requires `texlive-luatex` (fontspec + luaotfload) and one of the fonts above.

## Building a panel

```
cd figN && python figNx_<name>.py
```

Every panel script imports `nm_style.py` (fonts, type scale, palette, export settings)
and, for data panels, `remote_data.py` or the figure's own `figN_data.py`. Each writes
a vector PDF and a 600 dpi PNG preview. Extended Data Fig. 1 is built by
`../../scripts/figures/draw_ed_fig1.py`, which writes all three of its output copies.

## Style contract

- One sans family per figure; `nm_style.resolve_sans()` raises rather than falling back
  to DejaVu Sans, which is neither Arial- nor Helvetica-metric.
- House gene palette (`nm_style.GENE_COLORS`) means "gene" and is used only where the
  variable really is the gene. External atlases take the neutral slate ramp
  (`nm_style.DATASET_COLORS`), shared by Fig. 6 and Extended Data Fig. 1.
- One visual encoding per variable per figure; one key per figure at most, otherwise
  direct labels.
- Composites contain zero embedded raster: `pdfimages -list figN_composite.pdf | tail -n +3 | wc -l` -> 0.

## Verifying a rebuild

```
pdffonts   latex/figures/figN.pdf   # one family only
pdfinfo    latex/figures/figN.pdf   # page size unchanged
pdfimages -list latex/figures/figN.pdf | tail -n +3 | wc -l   # 0
```

For a change that is meant to be visual-only or style-only, render the old and new PDFs
at 300 dpi and diff the pixels; the differing regions should fall exactly where the
change was intended and nowhere else.

## Superseded material

`_superseded/` holds panels replaced during the Nature Methods compliance passes,
including externally authored art that was placed at scales of 0.18-0.41 with type
collapsing to 3.4 pt. It is kept for provenance only; do not reuse it, and do not
apply the current style contract to it.
