# Figure 1 — AllelePerturb defines protein-coding variant-level perturbation prediction

**Role:** Task definition + benchmark overview + evaluation protocol. No results (results begin at Fig 2).

## Files
| File | Type | Use |
|------|------|-----|
| `fig1_final.png` | GPT illustration (from design spec) | **Main manuscript figure** |
| `fig1_benchmark_overview.png` | Data-driven (matplotlib) | Draft / data-panel reference (panels a-c) |
| `fig1_benchmark_overview.pdf` | Data-driven vector | LaTeX-embeddable data version |

## Panels (in `fig1_final.png`)
- **a** Gene-level averaging hides allele-specific trajectories (schematic)
- **b** Lollipop: 472 variants across four protein architectures (TP53/KRAS/GATA1/JAK1)
- **c** Per-variant depth varies widely (98/93/255/26 variants; 52.9k/61.3k/149.2k/5.6k cells)
- **d** theta encodes protein-variant properties beyond gene identity (6-dim + PCA)
- **e** Evaluation protocol: direction recovery vs allele-specific ranking + 6 splits

## Production
- **Data version** (`fig1_benchmark_overview.png`): produced from `data/allele_perturb_bench.csv`.
  Design + draft is data-direct; final artwork rendered by GPT from the Claude design spec.
- Figure 1 is a **design-spec + GPT-art** figure (not pure data-direct). The data panels (b lollipop,
  c depth) are backed by the benchmark CSV; the schematic panels (a, d, e) are illustration.

## Data dependencies
- `data/allele_perturb_bench.csv` (472 variants, θ features, n_cells, splits)
