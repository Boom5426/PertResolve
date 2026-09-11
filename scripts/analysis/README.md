# Analysis generators

The scripts in this directory generate the canonical benchmark and
measurement-resolution tables. Committed outputs live under
`results/canonical/` and are kept separate from scratch reruns.

These are paper-specific generators, not alternate defaults for the generic
`pertresolve.resolution` API. In particular, the paper's six model heads across
three input representations, train-only standardization, same-gene candidate
pools and disjoint sampling protocol are preserved here. The generic API's
alpha-like predictor family and four-group diagnostic are useful for software
checks, but its `Δalpha`-like resolution is not the paper's observed `ΔPDS`.
Split-half outputs are described as empirical reproducibility references; they
are not hard ceilings or bounds and do not, by themselves, certify model-ranking
resolution.

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
