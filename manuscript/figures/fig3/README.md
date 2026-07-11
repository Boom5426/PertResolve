# Figure 3 — assembly guide

**One-line message:** split-half analysis reveals an allele-resolution measurement
floor (within-variant noise vs variant signal); TP53/KRAS/GATA1 sit near the floor,
JAK1 retains a wide window. Source of truth: `manuscript/latex/AllelePerturb_manuscript.tex`.

Palette: house set (TP53 `#5185C0`, KRAS `#E99D4E`, GATA1 `#8281B9`, JAK1 `#55966B`)
via the shared `manuscript/figures/nm_style.py`.

## Panel status (this pass)

| Panel | Content | Status | File / source |
|-------|---------|--------|---------------|
| a | split-half window schematic | **done** (programmatic vector) | `fig3a_schematic.py` |
| c | per-variant D_self/D_null (hero) | **done** | `fig3c_ratio.py` — `second_probe_rankability_table.csv` (native-depth) + `bootstrap_CIs.json` |
| d | effect-size × depth landscape | **done** | `fig3d_landscape.py` — native rows (effect_size, D_null, n_cells) |
| e | depth titration | **done** | `fig3e_titration.py` — `split_half_power_curve.csv` (PCA-50) |
| f | rankable fraction (4 allele genes) | **done** | `fig3f_rankable.py` — `unrankable_canonical.json` |
| b | representative closed/open distributions | **data done** | `fig3b_windows.py` — `_remote/unified/fig3b_selfnull_dist.csv` (200-seed split-half from grid_cbv on the server: TP53 P222P closed R=0.93, JAK1 R108Q open R=0.16) |
| g | detection vs identification (sibling-pair resolvable 0/0/3.8/78.5%) | **HELD** | no committed local file (remote) |
| h | model-free classifier AUROC (0.49/0.50/0.52/0.87) | **HELD** | no committed local file (remote) |

Run any panel: `python fig3<x>_*.py` (all import the shared `nm_style` and, for data
panels, `fig3_data.py`). Each writes a vector PDF (0 embedded raster) + 600 dpi PNG preview.

## Numbers (all verified to trace to committed files)

- D_self/D_null (mean, 95% CI): TP53 0.96 (0.95–0.98), KRAS 1.00 (0.99–1.02),
  GATA1 0.88 (0.86–0.90), JAK1 0.21 (0.12–0.33). Native-depth aggregation reproduces
  `canonical_numbers.json` exactly.
- Rankable at native depth: TP53 0%, KRAS 0%, GATA1 2.4%, JAK1 90%
  (n = 98 / 92 / 254 / 20 evaluated perturbations).
- Depth titration (PCA-50, fraction detectable): JAK1 = 1.00 at all rungs; TP53
  0.51→0.69 (no 300 rung), KRAS 0.57→0.74, GATA1 ~0.51–0.57 (plateau < 75%).
- Effect size (median): TP53 1.41, KRAS 1.11, GATA1 2.79, JAK1 13.04.

## Notes / decisions

- **External gene-level atlases moved out of Fig 3.** Per the agreed plan, Replogle
  (55.3%), Norman (3.4%), Adamson (14.6%) un-rankable belong in Fig 6 / Extended Data,
  not here; Fig 3 stays allele-focused. `unrankable_canonical.json` still holds those
  values for the ED/Fig 6 panel. The manuscript currently shows them in Fig 3d — that
  text/structure edit is pending sign-off.
- **3b/3g/3h are the strongest "ground-truth limitation" panels but have no committed
  local data** (computed remotely). Pull priority when the server returns: sibling-pair
  resolvability + single-cell classifier AUROC + representative per-cell distributions.
- 3a is programmatic (not AI) to stay vector and palette-consistent; swap for an AI
  schematic only if you prefer.
