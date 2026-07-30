#!/usr/bin/env python3
"""Allele-level datasets on the SAME resolution axis as the public benchmarks.

Same M=50 split-half, same fine alpha grid, same min_resolvable_gap coefficient as
benchmark_resolution.py, but loading the four allele genes from the grid arrays via
harness. Lets TP53/KRAS/GATA1/JAK1 be plotted on one axis with the 5 public datasets.
Isolated output -> <out>/allele_<gene>.json (does not touch manuscript).

Usage:
    python allele_resolution.py [--base BASE] --out OUT

``--base`` (or the ``VCCOMPASS_BASE`` environment variable) points at the compute
workspace holding ``unified/harness.py`` and ``allele_perturb_bench.csv``.
``--out`` is required and must not point inside this repository's ``results/``,
so a re-run can never overwrite a committed canonical table.
"""
import sys, json, os, numpy as np
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from alleleperturb.paths import resolve_base, add_harness_to_path, require_inputs

_parser = argparse.ArgumentParser(
    description="Allele-level resolution curves on the public-benchmark axis.")
_parser.add_argument("--base", default=None,
                     help="VCCompass compute workspace (default: $VCCOMPASS_BASE).")
_parser.add_argument("--out", required=True, type=Path,
                     help="Directory receiving the allele_<gene>.json files.")
_args = _parser.parse_args()

BASE = resolve_base(_args.base)
BENCH_CSV = BASE / "allele_perturb_bench.csv"
require_inputs(BENCH_CSV)
add_harness_to_path(BASE)
OUT_DIR = _args.out

import harness as H

WT = ('WT', 'wt', 'WT_control')
M = 50; NSEED = 8; NBOOT = 500
ALPHAS = [0.0, 0.5, 0.7, 0.85, 0.925, 0.96, 1.0]; NA = len(ALPHAS)
np.random.seed(0)


def energy(A, B):
    dab = np.sqrt(((A[:, None] - B[None]) ** 2).sum(-1)).mean()
    daa = np.sqrt(((A[:, None] - A[None]) ** 2).sum(-1)).mean()
    dbb = np.sqrt(((B[:, None] - B[None]) ** 2).sum(-1)).mean()
    return 2 * dab - daa - dbb


def pds_row(drow, ci, n):
    td = drow[ci]
    less = int((drow < td - 1e-12).sum()); eq = int((np.abs(drow - td) <= 1e-12).sum())
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1) if n > 1 else 0.5


def norm(Mx):
    return Mx / (np.linalg.norm(Mx, axis=1, keepdims=True) + 1e-12)


import pandas as pd
_bench = pd.read_csv(BENCH_CSV)
BENCH = {g: set(_bench[_bench.gene == g]["variant"]) for g in H.GENES}   # curated benchmark set (consistency)


def run(g):
    X, lab = H.load_gene(g); lab = np.asarray(lab)
    from sklearn.decomposition import PCA
    Xp = PCA(n_components=50, random_state=0).fit_transform(X.astype(np.float32))
    wt = np.where(np.isin(lab, WT))[0]
    perts = [v for v in np.unique(lab) if v not in WT and v in BENCH[g] and (lab == v).sum() >= 2 * M]
    cell = {p: Xp[lab == p] for p in perts}
    dself, dnull = {p: [] for p in perts}, {p: [] for p in perts}
    pv = {a: {p: [] for p in perts} for a in ALPHAS}
    for seed in range(NSEED):
        rng = np.random.RandomState(seed)
        be = Xp[wt[rng.choice(len(wt), M, replace=False)]].mean(0)
        bb = Xp[wt[rng.choice(len(wt), M, replace=False)]].mean(0)
        build, evald = {}, {}
        for p in perts:
            c = cell[p]; ii = rng.permutation(len(c)); A, B = c[ii[:M]], c[ii[M:2 * M]]
            build[p] = A.mean(0) - bb; evald[p] = B.mean(0) - be
            dself[p].append(energy(A, B))
            dnull[p].append(energy(A, Xp[wt[rng.choice(len(wt), M, replace=False)]]))
        cv = list(build.keys()); gm = np.stack([build[p] for p in cv]).mean(0)
        Tn = norm(np.stack([evald[p] for p in cv])); n = len(cv)
        for a in ALPHAS:
            Pred = norm(np.stack([(1 - a) * gm + a * build[p] for p in cv]))
            D = 1 - Pred @ Tn.T
            for i, p in enumerate(cv):
                pv[a][p].append(pds_row(D[i], i, n))
    cv = [p for p in perts if pv[1.0][p]]
    rank = [(np.median(dnull[p]) - np.median(dself[p])) > (np.percentile(dself[p], 97.5) - np.percentile(dself[p], 2.5)) for p in cv]
    pvm = {a: np.array([np.mean(pv[a][p]) for p in cv]) for a in ALPHAS}
    rng = np.random.RandomState(0); n = len(cv)
    gaps = [round(ALPHAS[k + 1] - ALPHAS[k], 3) for k in range(NA - 1)]; winhi = np.zeros(NA - 1); nc = 0
    for _ in range(NBOOT):
        bi = rng.randint(0, n, n); yy = [pvm[a][bi].mean() for a in ALPHAS]
        if all(yy[k + 1] > yy[k] for k in range(NA - 1)): nc += 1
        for k in range(NA - 1):
            if yy[k + 1] > yy[k]: winhi[k] += 1
    winhi /= NBOOT
    rg = [gaps[k] for k in range(NA - 1) if winhi[k] > 0.9]
    out = dict(dataset=f"allele_{g}", n_pert=n, frac_rankable=round(float(np.mean(rank)), 3),
               reliability_P_recover_order=round(nc / NBOOT, 3), oracle_ceiling=round(float(pvm[1.0].mean()), 3),
               min_resolvable_gap=(min(rg) if rg else None),
               gap_resolution={str(gaps[k]): round(float(winhi[k]), 3) for k in range(NA - 1)})
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(OUT_DIR / f"allele_{g}.json", "w"), indent=2)
    print(f"allele_{g:6s} n={n:3d} rankable={out['frac_rankable']:.2f} oracle={out['oracle_ceiling']:.2f} "
          f"min_gap={out['min_resolvable_gap']} gap_res={out['gap_resolution']}")


for g in H.GENES:
    run(g)
