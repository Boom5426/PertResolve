# Figure 4 — Allele-ranking failure persists across evaluation choices and model interfaces

**Role:** Robustness. The failure is invariant to distance metric, feature representation, and split.

## Files
| File | Use |
|------|-----|
| `fig4_composite.png` | **Main manuscript figure** (external-θ version) |
| `fig4_composite.pdf` | LaTeX-embeddable |

## Panels (built directly by draw script — no separate panel PNGs)
- **a** PDS across 5 splits (heatmap): near/below chance everywhere
- **b** Pearson-δ across splits: consistently positive
- **c** Ranking failure invariant to distance function (PDS cosine 0.468, L1 0.454, L2 0.454)
- **d** Dissociation holds across all genes and splits
- **e** Richer features do not close the gap (θ 0.45, ESM 0.45, ESM+θ 0.46)
- **f** Existing methods assume gene-level identity (CPA/CellFlow/STATE/GEARS/scDFM/Biolord failure modes)

## Production
- **Script:** `scripts/figures/draw_fig4.py`
- **Run:** `cd scripts/figures && python draw_fig4.py`
- **Loader:** `load_v4()` → `results/results_v4_10metrics.csv` (canonical: exttheta)

## Data dependencies
- `results/results_v4_exttheta.csv` (external-θ, de-leaked — canonical)
- `results/results_v4_10metrics.csv`
