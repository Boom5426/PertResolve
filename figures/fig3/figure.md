# Figure 3 — Split-half analysis reveals a narrow allele-level measurement window

**Role:** Mechanism. WHY ranking fails — measurement floor / split-half noise, not model weakness.

## Files
| File | Use |
|------|-----|
| `fig3_composite.png` | **Main manuscript figure** (CI version) |
| `fig3_composite.pdf` | LaTeX-embeddable |
| `panels/fig3b-e.png` | Individual data panels |

## Panels
- **a** Concept: D_self / D_null / signal window (schematic — GPT/placeholder, not data-direct)
- **b** Replicate noise approaches variant signal: D_self/D_null = TP53 0.96, KRAS 1.01, GATA1 0.90, JAK1 0.14
- **c** Depth × effect size jointly set the window (cells/variant vs ratio)
- **d** Detection floor extends to gene-level benchmarks (Replogle 55%, Norman 3%, Adamson 15%) — **DIAGNOSTIC**
- **e** Detection rate rises with depth, gene-dependently (n=50-300 subsample)

## Production
- **Script:** `scripts/figures/draw_fig3.py`
- **Run:** `cd scripts/figures && python draw_fig3.py`
- **Loaders:** `load_rankability()` → `second_probe_rankability_table.csv`; `load_power()` → `split_half_power_curve.csv`
- **Fig 3a concept schematic** needs GPT redraw (placeholder).

## Data dependencies
- `results/second_probe_rankability_table.csv` (55,548 rows; native + subsampled edist/pca/mmd)
- `results/split_half_power_curve.csv` (detection rate vs depth)
- `results/unrankable_canonical.json`, `results/canonical_numbers.json`

## Note vs Fig 5d
Fig 3d = **diagnostic** ("is the floor allele-specific or does it extend to gene-level?"). Fig 5d = **prescriptive** (depth-normalization comparison). Cross-referenced in manuscript.
