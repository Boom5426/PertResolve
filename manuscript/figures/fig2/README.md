# Figure 2 — assembly guide

**One-line message:** variant-level prediction reveals a direction-ranking dissociation —
every predictor recovers the perturbation direction but none ranks the correct allele.
Source of truth: `manuscript/latex/AllelePerturb_manuscript.tex`. House palette + shared
`nm_style.py`; data via `remote_data.py`.

## Panel status

| Panel | Content | Type | File / source |
|-------|---------|------|---------------|
| a | unified-evaluation schematic | AI schematic | `fig2a_prompt.md` |
| b | PDS forest (20 in-house predictors) | **data done** | `fig2b_pds_forest.py` — `_remote/unified/definitive_summary.csv` |
| c | Pearson-δ forest | **data done** | `fig2c_pearson_forest.py` — `results/results_v4_exttheta.csv` |
| d | PDS-vs-Pearson scatter | **data done** | `fig2d_scatter.py` — PDS from definitive_summary, Pearson from exttheta |
| e | DE fidelity gradient | **data done** | `fig2e_de_gradient.py` — exttheta; drawn with reproducible 65/0.25/0.15, manuscript text updated |
| f | per-gene dissociation gap | **data done** | `fig2f_gene_gap.py` — exttheta |
| g | representative TP53 variants | **data done** | `fig2g_pervariant.py` — exttheta (Ridge-esm TP53) |

Feature-space colour (2b/2c/2d) uses the shared `S.FEATURE_COLORS` (theta=blue, ESM=orange,
ESM+theta=green, reference=grey) — one consistent mapping across the three panels.

## Verified numbers

- **2b** PDS (multi-seed, `definitive_summary.csv`): 18 heads + Gene-mean + WT-null, range
  0.487-0.517, every 95% CI crosses 0.50 (external SOTA excluded -> Fig 4g).
- **2c** Pearson-δ (`results_v4_exttheta.csv`): per-method mean 0.555-0.647, all CIs > 0.
  NOTE: this is slightly below the manuscript's stated 0.60-0.68 (the three KNN heads sit at
  0.555-0.567); reported honestly. Text 0.60-0.68 may want widening to ~0.55-0.65.
- **2d**: all predictors cluster in the "direction without ranking" region (PDS ~0.5, Pearson ~0.6); WT-null at (0.50, 0).
- **2f** gap = mean(Pearson) - mean(PDS): TP53 0.31, KRAS 0.19, GATA1 0.13, JAK1 -0.04 (~ manuscript 0.30/0.18/0.15/0.02; JAK1 ~ 0).
- **2g** Ridge-esm TP53: Pearson stable 0.72-0.84, PDS 0.00-0.98 across held-out variants.

## RESOLVED: panel 2e (DE gradient)

The manuscript DE gradient 78/46/28 reproduced from **no committed or remote file**
(results_v4 gives ~65/0.25/0.15; no method/aggregation on the server reaches 78/46/28).
Per the author's decision (option A), 2e is drawn with the reproducible in-house exttheta
values **direction 65% / DE-LFC Spearman 0.25 / DE overlap 15%**, and the manuscript text
was corrected in two places (Results + Fig 2 caption). Same decreasing gradient, lower
absolute values.
