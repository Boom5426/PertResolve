# Project Structure — AllelePerturb

A reproducible map of the AllelePerturb benchmark: every figure ties back to its
source script, its input data, and the manuscript. Organized figure-by-figure for reuse.

```
AllelePerturb/
├── alleleperturb/            # pip-installable package (metrics, PDS, eval)
│   ├── bench.py              # benchmark loader
│   ├── features.py           # θ feature construction
│   ├── metrics.py            # top-level metric API
│   └── evaluation/           # pds.py, de_metrics.py, direction_metrics.py
├── data/                     # benchmark tables (inputs)
│   ├── allele_perturb_bench.csv          # 472 variants, θ, n_cells, splits (primary)
│   ├── allele_perturb_bench_v2.csv       # corrected external-hotspot version
│   ├── allele_perturb_bench_extθ.csv     # external-θ (de-leaked)
│   └── hotspot_external_definition.txt   # external hotspot criteria (COSMIC/IARC/ClinVar)
├── results/                  # computed result tables (see results/README.md)
├── scripts/figures/          # figure-generating code (run from this dir)
│   ├── fig_config.py         # shared style + data loaders (anchors DATA_DIR=results/)
│   ├── draw_fig2.py … draw_fig5.py
│   ├── all_splits_v4.py      # the full 20-method × 5-split × 4-gene grid
│   └── rankability_audit.py  # split-half D_self/D_null + energy distance
├── figures/                  # FIGURE-CENTRIC organization (one dir per figure)
│   ├── fig1/  … fig5/         # each: composite (.png+.pdf), panels/, figure.md manifest
│   └── extended_data/         # ED Fig 1 rankability sensitivity
└── manuscript/
    ├── AllelePerturb_manuscript_en.md   # canonical English manuscript (with equations)
    ├── AllelePerturb_manuscript_cn.md   # Chinese draft (superseded; EN is canonical)
    └── latex/
        ├── AllelePerturb_manuscript.tex # NM-style LaTeX (Palatino, lineno, natbib)
        ├── AllelePerturb_manuscript.pdf # compiled, 17 pp, all 5 figures embedded
        └── figures/          # fig1.pdf … fig5.pdf (embedded by \includegraphics)
```

## Figure → Script → Data map

| Figure | Role | Script | Input data | Type |
|--------|------|--------|-----------|------|
| **Fig 1** | Task + benchmark + protocol | (design spec) | `data/allele_perturb_bench.csv` | GPT art (data panels backed) |
| **Fig 2** | Direction-ranking dissociation | `draw_fig2.py` | `results/results_v4_10metrics.csv` + `bootstrap_CIs.json` | Data-direct |
| **Fig 3** | Measurement-window mechanism | `draw_fig3.py` | `second_probe_rankability_table.csv` + `split_half_power_curve.csv` | Data-direct (3a schematic) |
| **Fig 4** | Robustness (metric/feature/split) | `draw_fig4.py` | `results_v4_exttheta.csv` | Data-direct |
| **Fig 5** | Power-aware reporting protocol | `draw_fig5.py` | `second_probe_predictor_results.csv` + rankability table | Data-direct (5e workflow) |
| **ED Fig 1** | Rankability criterion sensitivity | (round-2 probe) | `results/rankability_sensitivity.csv` | Data-direct |

## Reproducing figures

```bash
cd scripts/figures
python draw_fig2.py     # -> ../../figures/composites/... (see --out flag)
python draw_fig3.py
python draw_fig4.py
python draw_fig5.py
```

Loaders in `fig_config.py` resolve `results/` automatically (anchored to `DATA_DIR`).
Each `figures/figN/figure.md` documents that figure's panels, script, and data dependencies.

## Compiling the manuscript

```bash
cd manuscript/latex
pdflatex AllelePerturb_manuscript.tex   # run twice for cross-references
```
Requires Palatino fonts (`tlmgr install mathpazo palatino psnfss` on TinyTeX).
Figure PDFs live in `manuscript/latex/figures/`.

## Canonical numbers

Frozen in `results/canonical_numbers.json`. Key values:
- 472 variants, 321,043 cells, 4 genes, 3 technologies
- D_self/D_null: TP53 0.96, KRAS 1.01, GATA1 0.90, JAK1 0.14
- Un-rankable (native): TP53/KRAS/GATA1/JAK1 = 100/100/98/10%
- No method exceeds chance PDS (0.5); Pearson-δ 0.60–0.69
