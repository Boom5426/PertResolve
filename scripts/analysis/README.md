# Analysis generators

The scripts in this directory generate the canonical benchmark and
measurement-resolution tables. Committed outputs live under
`results/canonical/` and are kept separate from scratch reruns.

Representative generators include:

| Script | Main output or role |
|---|---|
| `train_only_grid.py` | feature-based model evaluation and PDS summaries |
| `oracle_ceiling.py` | split-half reproducibility reference |
| `resolution_panel_v2.py` | the broader PertResolve-Bench panel |
| `resolution_scaling.py` | controlled resolution and depth analyses |
| `resolved_pair_accuracy.py` | model discrimination on measurement-resolved pairs |
| `fig2_topk.py` and `split_half_topk.py` | top-k recovery summaries |
| `harness.py` | shared pseudobulk targets, splits and scoring conventions |

## Running analyses

Analyses that need single-cell matrices read a processed-data workspace that
is not redistributed with the repository. Name it explicitly; it is never
guessed:

```bash
export PERTRESOLVE_DATA=/path/to/processed-data  # or pass --base
python scripts/analysis/oracle_ceiling.py --out /scratch/pertresolve-rerun
```

The output path is required and is refused if it resolves inside this
repository's `results/` directory. Compare scratch outputs with the committed
tables before promoting anything into the public results tree.
