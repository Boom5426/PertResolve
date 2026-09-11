# Reviewer controls

Robustness and sanity controls that sit alongside the main evaluation grid. Each table is
small and committed so the Supplementary Note values regenerate from a tracked file.

| File | Backs | Script |
|---|---|---|
| `rep_grid_summary.csv` | Supplementary Note "Representation robustness": PDS-cosine and Pearson-delta per variant representation, single deterministic pseudobulk draw | `rep_grid.py` |
| `am_summary.csv` | the AlphaMissense rows of the same range (`alphamissense`, `alphamissense+theta`) | `am_analysis.py` |
| `am_effect_spearman.csv` | AlphaMissense pathogenicity versus measured effect size, per gene | `am_analysis.py` |
| `residual_decomp_measurement.csv` | gene-shared versus allele-specific decomposition, per gene, with full and residual oracle | `residual_decomp.py` |
| `residual_decomp_models.csv` | full and residual Pearson-delta per predictor | `residual_decomp.py` |
| `oracle_sensitivity.csv` | oracle ceiling at 300 cells, at full depth, and for a perfect prediction scored against noisy truth | `oracle_sensitivity.py` |
| `residue_leakage.csv` | training-to-held-out residue overlap per split (47 / 0 / 10 / 37 / 54 per cent) | `residue_independence.py` |
| `residue_independence_summary.csv` | residue-cluster bootstrap per predictor, random split | `residue_independence.py` |
| `wt_control_diagnostic.csv` | wild-type-versus-wild-type control for the split-half window | (remote metrology harness) |
| `esm2_extract.py` | ESM2 embedding extraction for the four ESM2 constructions | |

## Running the scripts

`rep_grid.py` and `am_analysis.py` need the large arrays in the compute workspace, so they take
the workspace root as an argument rather than hard-coding it:

```bash
python rep_grid.py    --base /path/to/processed-data --out /scratch/rep-grid
python am_analysis.py --base /path/to/processed-data --out /scratch/am-analysis
# or: export PERTRESOLVE_DATA=/path/to/processed-data
```

## Two caveats that must not be lost

**1. The MLP head is constructed by keyword, and this is load-bearing.** In scikit-learn 1.7
and later the first positional parameter of `MLPRegressor` is `loss`, not `hidden_layer_sizes`.
The original `MLPRegressor((128, 64), ...)` therefore raised `InvalidParameterError` on every
fit, and a bare `except Exception: continue` discarded it, so the representation grid silently
ran on five heads instead of six. `rep_grid_summary.csv` now carries `n_heads_scored`, `heads`
and `n_head_failures`, and both scripts count head-fit failures per gene, split and head, so
the same silent omission cannot recur. The superseded five-head table and the full before and
after comparison are in `../deprecated/rep_grid_summary_5heads.csv` and its README entry.

**2. Do not report the pooled AlphaMissense correlation.** `am_effect_spearman.csv` includes a
`POOLED` row (-0.3798) that disagrees in sign and magnitude with every per-gene value (TP53
-0.0752, KRAS -0.0538, GATA1 0.0874, JAK1 0.4091). Cross-gene pooling is not interpretable
here, because both the AlphaMissense pathogenicity distribution and the scale of the measured
transcriptional effect differ between genes, so the pooled statistic is driven by between-gene
offsets rather than by any within-gene relationship. The row is kept for provenance only. The
JAK1 value should also not be read as a positive result: only 11 of its 26 variants carry an
AlphaMissense score. The manuscript reports the within-gene TP53, KRAS and GATA1 values only.

## Variant universes

`unified/real_deltas.npz` in the compute workspace holds 97 KRAS keys, five more
(`AG11TD`, `AG59GV`, `C185Y`, `K179R`, `M170L`) than the curated benchmark's 92. Those five are
absent from `pertresolve_bench.csv`, which is what `harness.split_vars` reads, so they never
enter a train, test or candidate set, and none of them carries an AlphaMissense score. Every
number in these tables is on the benchmark universe of 470 variants (TP53 98, KRAS 92,
GATA1 254, JAK1 26).
