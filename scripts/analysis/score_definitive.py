#!/usr/bin/env python3
"""Canonical multi-seed scoring of every method in the PertResolve comparison.

This is the scoring path behind the main-text model comparison. Each held-out
variant's PDS is averaged over ``NSEED`` independent pseudobulk subsamples before
any pooling, which removes the arbitrariness of a single stochastic draw. The
point estimate is then the mean over distinct variants, each variant's score
first averaged over the splits that hold it out; the cell-wise pooling described
below is reported alongside it rather than as the estimate.

Two resampling units are reported, because they answer different questions and
the difference between them is not negligible:

``variant`` (primary)
    the bootstrap resamples the distinct (gene, variant) pairs. A variant held
    out by three splits contributes one entry, not three, so the interval
    reflects the number of alleles actually measured. This is the unit the
    manuscript's Study design declares.

``cell`` (secondary)
    the bootstrap resamples held-out variants within each (split, gene) cell,
    which weights a variant by how many splits hold it out. Reported for
    comparability with the per-split grid, and because it is the convention the
    earlier runs of this study used.

The per-seed, per-variant scores are written out as well, so that the pooling
above can be recomputed, and so that subsets of the benchmark (for instance the
single-substitution missense conditions alone) can be scored without rerunning
the grid.
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import reject_repo_results, require_inputs, resolve_base

_ap = argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--base", default=None,
                 help="directory holding the gene arrays and preds5/ "
                      "(env: PERTRESOLVE_DATA)")
_ap.add_argument("--preds", default=None,
                 help="directory of per-method prediction .npz files "
                      "(default: <base>/unified/preds5)")
_ap.add_argument("--out", required=True, metavar="OUT_DIR",
                 help="directory receiving the tables; may not be inside results/")
_ap.add_argument("--tag", default="",
                 help="suffix distinguishing a variant run, e.g. 'trainonly'")
_args = _ap.parse_args()

BASE = resolve_base(_args.base)
OUT_DIR = reject_repo_results(_args.out)
OUT_DIR.mkdir(parents=True, exist_ok=True)
PREDS = Path(_args.preds) if _args.preds else BASE / "unified" / "preds5"
require_inputs(PREDS)
SUFFIX = f"_{_args.tag}" if _args.tag else ""

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness as H  # noqa: E402  (import after sys.path is prepared)

H.set_base(BASE)

np.random.seed(0)
NSEED = 15
NSUB = 300
NBOOT = 2000
SPLITS = ['split1', 'split2', 'split3', 'split5', 'split6']

cand = {(g, s): [v for v in (H.split_vars(g, s)[0] + H.split_vars(g, s)[1])]
        for g in H.GENES for s in SPLITS}
gene_cells = {g: (lambda X, lab: (X, np.asarray(lab)))(*H.load_gene(g)) for g in H.GENES}


def real_deltas_seed(seed):
    """Per-variant pseudobulk profiles for one subsample seed.

    Args:
        seed: index in ``range(NSEED)``; the per-gene ``RandomState`` is
            ``1000 + 7 * seed + gene_seed``, so genes stay independent and the
            draw is reproducible.

    Returns:
        ``{"<gene>__<variant>": delta}`` over variants carrying at least five cells.
    """
    rd = {}
    for g in H.GENES:
        X, lab = gene_cells[g]
        rng = np.random.RandomState(1000 + 7 * seed + H._GENE_SEED[g])
        wt = np.where(np.isin(lab, H.WT_TAGS))[0]
        wt = rng.choice(wt, NSUB, replace=False) if len(wt) > NSUB else wt
        wm = X[wt].mean(0)
        for v in np.unique(lab):
            if v in H.WT_TAGS:
                continue
            idx = np.where(lab == v)[0]
            if len(idx) < 5:
                continue
            if len(idx) > NSUB:
                idx = rng.choice(idx, NSUB, replace=False)
            rd[f"{g}__{v}"] = (X[idx].mean(0) - wm).astype(np.float32)
    return rd


def norm(M):
    """Row-normalize, leaving a zero row at zero rather than dividing by zero."""
    return M / (np.linalg.norm(M, axis=1, keepdims=True) + 1e-12)


def pds_row(drow, ci, n):
    """Tie-aware mid-rank PDS for one prediction against ``n`` candidates.

    A non-finite distance to the true match means the prediction could not be
    ranked at all; it takes the chance value, which is the harness convention
    the manuscript reports alongside the exclusion alternative.
    """
    td = drow[ci]
    if not np.isfinite(td):
        return 0.5
    less = int((drow < td - 1e-12).sum())
    eq = int((np.abs(drow - td) <= 1e-12).sum())
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1) if n > 1 else 0.5


methods = [os.path.basename(f)[:-4] for f in sorted(glob.glob(str(PREDS / "*.npz")))]
if not methods:
    raise SystemExit(f"No prediction files found under {PREDS}")
P = {m: np.load(PREDS / f"{m}.npz") for m in methods}

# acc[method][(gene, split)][variant] = list of per-seed PDS
acc = {m: {} for m in methods}
long_rows = []
for seed in range(NSEED):
    real = real_deltas_seed(seed)
    for m in methods:
        for g in H.GENES:
            for s in SPLITS:
                cv = [v for v in cand[(g, s)] if f"{g}__{v}" in real]
                idx = {v: i for i, v in enumerate(cv)}
                te = [v for v in H.split_vars(g, s)[1]
                      if f"{g}__{s}__{v}" in P[m].files and v in idx]
                if not te:
                    continue
                Pm = norm(np.stack([P[m][f"{g}__{s}__{v}"] for v in te]))
                Tm = norm(np.stack([real[f"{g}__{u}"] for u in cv]))
                D = 1 - Pm @ Tm.T
                for i, v in enumerate(te):
                    score = pds_row(D[i], idx[v], len(cv))
                    acc[m].setdefault((g, s), {}).setdefault(v, []).append(score)
                    long_rows.append((m, g, s, v, seed, score))
    print(f"seed {seed} done", flush=True)

long = pd.DataFrame(long_rows, columns=["method", "gene", "split", "variant", "seed", "pds"])
long_path = OUT_DIR / f"per_seed_variant_pds{SUFFIX}.csv.gz"
long.to_csv(long_path, index=False, float_format="%.6f", compression="gzip")
print(f"wrote {len(long):,} per-seed per-variant scores -> {long_path}", flush=True)


def bootstrap_cells(cells, rng):
    """Resample held-out variants within each (split, gene) cell.

    Args:
        cells: ``{(gene, split): {variant: seed_averaged_pds}}``.
        rng: the shared generator, so the runs are reproducible.

    Returns:
        Point estimate and the 2.5 / 97.5 percentiles over ``NBOOT`` resamples.
    """
    cellmeans = {k: np.mean(list(d.values())) for k, d in cells.items()}
    obs = float(np.mean(list(cellmeans.values())))
    boot = np.empty(NBOOT)
    arrays = [np.array(list(d.values())) for d in cells.values()]
    for b in range(NBOOT):
        boot[b] = np.mean([a[rng.integers(0, len(a), len(a))].mean() for a in arrays])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return obs, float(lo), float(hi)


def bootstrap_variants(cells, rng):
    """Resample distinct (gene, variant) pairs, averaging a variant's splits first.

    This is the clustered interval: repeated scorings of one allele under
    different splits are not independent, so they are collapsed before
    resampling rather than treated as separate draws.
    """
    per_var = {}
    for (g, _s), d in cells.items():
        for v, score in d.items():
            per_var.setdefault((g, v), []).append(score)
    vals = np.array([np.mean(x) for x in per_var.values()])
    obs = float(vals.mean())
    boot = np.array([vals[rng.integers(0, len(vals), len(vals))].mean()
                     for _ in range(NBOOT)])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return obs, float(lo), float(hi), len(vals)


rows = []
for m in methods:
    cells = {k: {v: float(np.mean(l)) for v, l in d.items()} for k, d in acc[m].items()}
    rng = np.random.default_rng(0)
    obs_c, lo_c, hi_c = bootstrap_cells(cells, rng)
    rng = np.random.default_rng(0)
    obs_v, lo_v, hi_v, n_var = bootstrap_variants(cells, rng)
    rows.append(dict(
        method=m,
        PDS_variant=round(obs_v, 4), variant_lo=round(lo_v, 4), variant_hi=round(hi_v, 4),
        n_variants=n_var, variant_crosses=bool(lo_v <= 0.5 <= hi_v),
        PDS_cell=round(obs_c, 4), cell_lo=round(lo_c, 4), cell_hi=round(hi_c, 4),
        n_cells=len(cells), cell_crosses=bool(lo_c <= 0.5 <= hi_c),
    ))

df = pd.DataFrame(rows).sort_values("PDS_variant", ascending=False)
out_path = OUT_DIR / f"definitive_summary{SUFFIX}.csv"
df.to_csv(out_path, index=False)
print(df.to_string(index=False))
print(f"\nwrote {out_path}")

heads = df[df.method.str.contains("Ridge|Lasso|RF|GBoost|KNN|MLP")]
print(f"\n18 feature-model heads: variant-unit mean {heads.PDS_variant.mean():.4f}, "
      f"cell-unit mean {heads.PDS_cell.mean():.4f}")
print(f"all methods variant-unit CI crosses 0.5: {df.variant_crosses.all()}; "
      f"exceptions: {list(df[~df.variant_crosses].method)}")
print(f"all methods cell-unit CI crosses 0.5: {df.cell_crosses.all()}; "
      f"exceptions: {list(df[~df.cell_crosses].method)}")
