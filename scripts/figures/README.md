# AllelePerturb Figure Drawing Scripts

## Setup

```bash
cd /home/boom/ICLR/AllelePerturb/scripts/figures
# Requires: numpy, pandas, matplotlib, scikit-learn
```

## Usage

```bash
# Draw all figures (uses data from ../../results/)
python draw_fig2.py
python draw_fig3.py
python draw_fig4.py
python draw_fig5.py

# Custom data path
python draw_fig2.py --data /path/to/results_v4_10metrics.csv

# Custom output
python draw_fig2.py --out /path/to/output.png
```

## Files

| File | Description |
|------|-------------|
| `fig_config.py` | Shared config: colors, styles, data loaders, helpers |
| `draw_fig2.py` | Figure 2: Direction-ranking dissociation (5 panels) |
| `draw_fig3.py` | Figure 3: Split-half measurement window (5 panels) |
| `draw_fig4.py` | Figure 4: Cross-split/metric/feature robustness (6 panels) |
| `draw_fig5.py` | Figure 5: Rankability prediction + design guidance (5 panels) |

## Data dependencies (in ../../results/)

| File | Rows | Source |
|------|------|--------|
| `results_v4_10metrics.csv` | 10,276 | 20 methods × 5 splits × 4 genes × 10 metrics |
| `second_probe_rankability_table.csv` | 55,548 | Split-half D_self/D_null per perturbation |
| `split_half_power_curve.csv` | 30 | Detection rate vs subsampled cells |
| `second_probe_predictor_results.csv` | 7 | LODO AUROC + calibration per dataset |

## Editing tips

- All colors defined in `fig_config.py` GENE_COLORS / FEAT_COLORS / BASELINE_COLOR
- Font sizes: title=7, axis=7, tick=6, legend=5.5, annotation=5
- Panel letters: `panel_letter(ax, 'a')` — bold, top-left
- To change a single panel: find the `# === PANEL X ===` block in draw_figN.py
