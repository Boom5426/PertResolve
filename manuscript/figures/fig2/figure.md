# Figure 2 — All methods recover perturbation direction but cannot rank variant identity

**Role:** The core phenomenon. Direction-ranking dissociation across 20 method x feature combinations.

## Files
| File | Use |
|------|-----|
| `fig2_composite.png` | **Main manuscript figure** (CI version, 95% bootstrap whiskers) |
| `fig2_composite.pdf` | LaTeX-embeddable |
| `panels/fig2a-e.png` | Individual panels |

## Panels
- **a** PDS (cosine) + Pearson-δ per method, 95% CI whiskers — no method exceeds chance (0.5)
- **b** Direction-without-ranking regime scatter (all methods cluster left of chance line)
- **c** DE fidelity: direction agreement 78%, DE-LFC Spearman 0.46, DE overlap 28%
- **d** Per-gene dissociation gap: TP53 0.30 > KRAS 0.18 > GATA1 0.15 > JAK1 0.02
- **e** Per-variant: stable direction, unstable ranking (synonymous vs hotspot)

## Production
- **Script:** `scripts/figures/draw_fig2.py`
- **Run:** `cd scripts/figures && python draw_fig2.py`  (uses `from fig_config import *`)
- **Loader:** `load_v4()` → `results/results_v4_10metrics.csv`

## Data dependencies
- `results/results_v4_10metrics.csv` (10,276 rows = 20 methods × 5 splits × 4 genes × 10 metrics)
- `results/bootstrap_CIs.json` (95% CI whiskers)
- Canonical (external-θ, de-leaked): `results/results_v4_exttheta.csv`
