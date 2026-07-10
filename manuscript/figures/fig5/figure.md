# Figure 5 — A power-aware reporting protocol for allele-resolution perturbation prediction

**Role:** Prescription. Rankability is predictable from pilot data; recommended workflow.

## Files
| File | Use |
|------|-----|
| `fig5_composite.png` | **Main manuscript figure** |
| `fig5_composite.pdf` | LaTeX-embeddable |

## Panels (built directly by draw script)
- **a** Per-perturbation effect size predicts rankability: LODO AUROC + bootstrap CI (mean 0.96 over 5 evaluable datasets; TP53/KRAS uniformly un-rankable, not evaluable)
- **b** Rankability set by effect-size regime, not cell count alone
- **c** Floor is intrinsic to distributional evaluation at finite depth (Replogle/Norman/Adamson native vs matched-n50), **PRESCRIPTIVE**
- **d** Recommended power-aware allele-resolution workflow (5 steps)

## Production
- **Script:** `scripts/figures/draw_fig5.py`
- **Run:** `cd scripts/figures && python draw_fig5.py`
- **Loaders:** honest predictor read from `results/rankability_predictor_honest.csv` (built by `rankability_predictor.py`); `load_rankability()` → `second_probe_rankability_table.csv`. (`load_predictor()` is deprecated; it pointed to the superseded row-level artifact.)
- Panel e workflow boxes may be GPT-polished from the data-direct layout.

## Data dependencies
- `results/rankability_predictor_honest.csv` (honest per-perturbation LODO predictor)
- `results/second_probe_rankability_table.csv`
- `results/canonical_numbers.json`

## Note vs Fig 3d
Fig 5c = **prescriptive** (depth-normalization). Fig 3d = **diagnostic** (floor is allele-specific). Cross-referenced.
