# Figures

Figure-centric organization: each figure has its own directory containing the
composite (PNG for viewing + PDF for LaTeX), individual panels, and a `figure.md`
manifest documenting panels, source script, and data dependencies.

| Dir | Figure | Composite | Panels | Manifest |
|-----|--------|-----------|--------|----------|
| `fig1/` | Task + benchmark overview | `fig1_final.png` (GPT), `fig1_benchmark_overview.{png,pdf}` (data) | in composite | `figure.md` |
| `fig2/` | Direction-ranking dissociation | `fig2_composite.{png,pdf}` | `panels/fig2a-e.png` | `figure.md` |
| `fig3/` | Measurement window | `fig3_composite.{png,pdf}` | `panels/fig3b-e.png` | `figure.md` |
| `fig4/` | Robustness | `fig4_composite.{png,pdf}` | (built by script) | `figure.md` |
| `fig5/` | Reporting protocol | `fig5_composite.{png,pdf}` | (built by script) | `figure.md` |
| `extended_data/` | ED Fig 1 sensitivity | `ED_fig1_rankability_sensitivity.{png,pdf}` | — | `figure.md` |

## Canonical versions

All composites here are the **final review-round-2 versions**:
- Fig 2, 3: CI versions (95% bootstrap whiskers)
- Fig 4: external-θ (de-leaked) version
- Fig 1: GPT art is the manuscript final (`fig1_final.png`); data version kept for reference

## Regenerating

Data-driven figures (2-5) are produced by `../scripts/figures/draw_figN.py`.
See each `figN/figure.md` for the exact script + data inputs.
Fig 1 and the Fig 3a schematic are GPT-rendered from Claude design specs.
