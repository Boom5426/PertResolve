# Figure 5 — A power-aware reporting protocol for allele-resolution perturbation prediction

**Role:** Prescription. Rankability is predictable from pilot data; recommended workflow.

## Files
| File | Use |
|------|-----|
| `fig5_composite.png` | **Main manuscript figure** |
| `fig5_composite.pdf` | LaTeX-embeddable |

## Panels (built directly by draw script)
- **a** Effect size predicts rankability: LODO ROC (mean AUROC core 0.974, all 0.965)
- **b** Discrimination strong, absolute calibration weak (R²=0.11, mean|err|=26pp)
- **c** Rankability set by effect-size regime, not cell count alone
- **d** Floor is intrinsic to distributional evaluation at finite depth (Replogle/Norman/Adamson native vs matched-n50) — **PRESCRIPTIVE**
- **e** Recommended power-aware allele-resolution workflow (5 steps)

## Production
- **Script:** `scripts/figures/draw_fig5.py`
- **Run:** `cd scripts/figures && python draw_fig5.py`
- **Loaders:** `load_predictor()` → `second_probe_predictor_results.csv`; `load_rankability()` → `second_probe_rankability_table.csv`
- Panel e workflow boxes may be GPT-polished from the data-direct layout.

## Data dependencies
- `results/second_probe_predictor_results.csv` (LODO predictor, 7 datasets)
- `results/second_probe_rankability_table.csv`
- `results/canonical_numbers.json`

## Note vs Fig 3d
Fig 5d = **prescriptive** (depth-normalization). Fig 3d = **diagnostic** (floor is allele-specific). Cross-referenced.
