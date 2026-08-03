# Canonical analysis generators

Every table in `results/canonical/` is cited by a figure, a Supplementary table or a
number in the text. Until now six of them had no committed generator: the scripts lived
only in the compute workspace under `unified/`, so a reader could see the number but not
the code that produced it. These are those scripts, rewritten onto this repository's path
contract and verified to reproduce their tables byte for byte.

| Script | Produces | Backs |
|---|---|---|
| `oracle_ceiling.py` | `oracle_ceiling.csv` | Fig. 5b, 5c replicate ceiling; SI Table 11 |
| `pairwise_resolvability.py` | `pairwise_resolvability.csv` | printed Fig. 4d sibling identification; SI Table 11 |
| `classifier_two_sample.py` | `classifier_two_sample.csv` | printed Fig. 4e classifier AUROC |
| `metric_floor.py` | `metric_family_floor.csv` | SI Table 4 metric-family invariance |
| `controlled_predictors.py` | `controlled_recovery.csv` | Fig. 5f, 5g benchmark-resolution transition |
| `floor_law_fit.py` | `floor_law.csv` | Methods, analytical scaling |

`subspace_test.py` was already here and follows the same contract.

## Running them

None of these can run from a checkout alone. They read the per-gene expression arrays and
the shared scorer, which are not redistributed (see the manuscript's Data availability), so
the workspace holding them is named explicitly and never guessed:

```sh
export VCCOMPASS_BASE=/path/to/VCCompass          # or pass --base
python scripts/analysis/oracle_ceiling.py --out /scratch/rerun
```

`--out` is required and is refused if it resolves inside this repository's `results/`.
Those tables back manuscript numbers, and a re-run that wrote straight back into that tree
could replace a cited value with the output of a changed script. Write to scratch, diff
against the committed table, then promote the file deliberately.

## Verification, 2026-08-03

Both the workspace originals and these rewritten copies were run in the remote `Agent`
environment against the same inputs, writing to scratch directories. All six outputs are
byte-identical to the committed tables in `results/canonical/`, so the rewrite changed no
number and the tables are reproducible from a fixed seed.

The rewrite touched exactly three lines per script: the hard-coded workspace root, the
hard-coded output path, and the message naming it. The numerical bodies are line-identical
to the workspace originals.

## Two properties of these analyses worth knowing before citing them

**They read `allele_perturb_bench.csv`, not `allele_perturb_bench_v2.csv`.** The shared
scorer builds its theta vectors from the v1 table, so the canonical in-house model grid is
computed on v1 theta. The two tables differ only in `is_hotspot`, and only on 182 of 472
rows. In v1 that column is not the indicator the Methods describes: it takes 83 distinct
values between -0.428 and 1.0, and 108 rows carry a value that is neither 0 nor 1, all of
them in TP53 (84 of 98) and KRAS (24 of 93), none in GATA1 or JAK1. v2 replaces it with the
external annotation and preserves the original as `is_hotspot_leaked_OLD`.

What those 108 values are is not settled by the tables themselves. They do **not** track the
measured per-variant effect size (Spearman -0.15 for TP53 and +0.16 for KRAS within those
rows), so they do not appear to carry response information; they do track the biophysical
change terms (up to 0.47 with `d_charge`, 0.46 with `d_vol`), which are data-independent.
The practical consequence is narrow: the sixth theta component of the canonical grid is not
the external annotation the Methods states, but it is not a readout of the outcome either,
and every score involved sits at chance regardless. Figure 1 and
`results/results_v4_exttheta.csv` use v2.

**Two analyses cap their work for speed, and the caps are real.** `metric_floor.py` uses at
most 40 variants per gene (`MAXVAR`), drawn by a fixed seed, and `classifier_two_sample.py`
scores 100 random sibling pairs per gene (`NPAIR`) with 30 permuted controls. Both are
deterministic, but neither is exhaustive, so their per-gene values are estimates over a
sample of variants rather than over all of them.
