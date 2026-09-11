#!/usr/bin/env python3
"""Top-3 #1: gene-shared vs allele-specific variance decomposition.

delta_v = delta_bar_gene + r_v   (gene-shared direction + allele-specific residual)

Part A (measurement side, split-half, multi-seed):
  - variance decomposition: fraction of per-variant effect that is gene-shared vs residual
  - FULL oracle PDS   (rank d_A(v) among {d_B(u)})            -> existing oracle
  - RESIDUAL oracle PDS (rank rA(v) among {rB(u)}, r = d - mean_u d)  -> is the ALLELE-SPECIFIC
    residual even measurable? Removing the shared gene direction isolates allele signal.
Part B (model side): FULL Pearson-delta vs RESIDUAL Pearson-delta per model
  (does the model recover allele-specific direction, or only the gene-shared program?)

Usage:
  python residual_decomp.py --out /path/to/output_dir [--base /path/to/processed-data]
  PERTRESOLVE_DATA=/path/to/processed-data python residual_decomp.py --out /path/to/output_dir
"""
import argparse
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
from pertresolve.paths import (add_harness_to_path, reject_repo_results,  # noqa: E402
                                 require_inputs, resolve_base)

# Arguments are parsed before the scientific stack is imported, so --help works
# without numpy/scipy/pandas installed and a bad --base fails before any load.
ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--base", default=None,
                help="directory holding the shared scorerharness.py, "
                     "unified/real_deltas.npz, unified/preds5/ and "
                     "pertresolve_bench.csv (env: PERTRESOLVE_DATA)")
ap.add_argument("--out", required=True,
                help="directory to write residual_decomp_measurement.csv and "
                     "residual_decomp_models.csv into")
args = ap.parse_args()

import numpy as np, pandas as pd, glob, os  # noqa: E402
from scipy.stats import pearsonr  # noqa: E402

BASE = resolve_base(args.base)
OUT_DIR = reject_repo_results(args.out)
BENCH_CSV = BASE / "pertresolve_bench.csv"
REAL_NPZ = BASE / "unified" / "real_deltas.npz"
PREDS_DIR = BASE / "unified" / "preds5"
require_inputs(BENCH_CSV, REAL_NPZ, PREDS_DIR)

add_harness_to_path(BASE)
import harness as H  # noqa: E402

NSUB = 300; NSEED = 15
SPLITS = ['split1', 'split2', 'split3', 'split5', 'split6']
WT_TAGS = ('WT', 'wt', 'WT_control')
np.random.seed(0)
gene_cells = {g: (lambda X, l: (X, np.asarray(l)))(*H.load_gene(g)) for g in H.GENES}
_bench = pd.read_csv(BENCH_CSV)
BENCH = {g: set(_bench[_bench.gene == g]["variant"]) for g in H.GENES}   # curated benchmark set (consistency)


def pds_row(drow, ci, n):
    td = drow[ci]
    if not np.isfinite(td):
        return 0.5
    less = int((drow < td - 1e-12).sum()); eq = int((np.abs(drow - td) <= 1e-12).sum())
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1) if n > 1 else 0.5


def norm(Mx):
    return Mx / (np.linalg.norm(Mx, axis=1, keepdims=True) + 1e-12)


# ---------- Part A: variance decomposition + full/residual oracle ----------
print("=== Part A: variance decomposition + full vs residual oracle (per gene) ===")
rowsA = []
for g in H.GENES:
    X, lab = gene_cells[g]
    fshare, full_or, resid_or = [], [], []
    for seed in range(NSEED):
        rng = np.random.RandomState(1000 + 7 * seed + H._GENE_SEED[g])
        wt = np.where(np.isin(lab, WT_TAGS))[0]; rng.shuffle(wt); h = len(wt) // 2
        wmA = X[wt[:h][:NSUB]].mean(0); wmB = X[wt[h:][:NSUB]].mean(0)
        dA, dB, vs = {}, {}, []
        for v in np.unique(lab):
            if v in WT_TAGS or v not in BENCH[g]:
                continue
            idx = np.where(lab == v)[0]
            if len(idx) < 5:
                continue
            rng.shuffle(idx); hh = len(idx) // 2
            dA[v] = X[idx[:hh][:NSUB]].mean(0) - wmA
            dB[v] = X[idx[hh:][:NSUB]].mean(0) - wmB
            vs.append(v)
        if len(vs) < 3:
            continue
        DA = np.stack([dA[v] for v in vs]); DB = np.stack([dB[v] for v in vs])
        gmB = DB.mean(0)                      # gene-shared direction (from B half)
        rA = DA - DA.mean(0); rB = DB - gmB   # allele-specific residuals
        # variance decomposition (from B half): total = shared + residual (exact, residuals zero-mean)
        total = float((DB ** 2).sum(1).mean()); shared = float((gmB ** 2).sum()); resid = float((rB ** 2).sum(1).mean())
        fshare.append(shared / (total + 1e-12))
        # full oracle: rank d_A(v) among {d_B(u)}
        n = len(vs); Df = 1 - norm(DA) @ norm(DB).T
        full_or.append(np.mean([pds_row(Df[i], i, n) for i in range(n)]))
        # residual oracle: rank rA(v) among {rB(u)}
        Dr = 1 - norm(rA) @ norm(rB).T
        resid_or.append(np.mean([pds_row(Dr[i], i, n) for i in range(n)]))
    rowsA.append(dict(gene=g, n_var=len(vs),
                      frac_gene_shared=round(float(np.mean(fshare)), 3),
                      frac_allele_residual=round(1 - float(np.mean(fshare)), 3),
                      oracle_full=round(float(np.mean(full_or)), 3),
                      oracle_residual=round(float(np.mean(resid_or)), 3)))
dfA = pd.DataFrame(rowsA)
print(dfA.to_string(index=False))
OUT_DIR.mkdir(parents=True, exist_ok=True)
dfA.to_csv(OUT_DIR / "residual_decomp_measurement.csv", index=False)

# ---------- Part B: model FULL vs RESIDUAL Pearson-delta ----------
print("\n=== Part B: model FULL vs RESIDUAL Pearson-delta (pooled over genes) ===")
real = dict(np.load(REAL_NPZ))
methods = [os.path.basename(f)[:-4]
           for f in sorted(glob.glob(os.path.join(glob.escape(str(PREDS_DIR)), "*.npz")))]
rowsB = []
for m in methods:
    P = np.load(PREDS_DIR / f"{m}.npz")
    full_p, resid_p = [], []
    for g in H.GENES:
        # collect held-out test variants (union over splits) with both truth and prediction
        tv, T, Pm = [], [], []
        for s in SPLITS:
            for v in H.split_vars(g, s)[1]:
                k = f"{g}__{v}"; pk = f"{g}__{s}__{v}"
                if k in real and pk in P.files and v not in tv:
                    tv.append(v); T.append(real[k]); Pm.append(P[pk])
        if len(tv) < 3:
            continue
        T = np.stack(T); Pm = np.stack(Pm)
        tgm = T.mean(0); pgm = Pm.mean(0)          # gene means (true / predicted)
        rT = T - tgm; rP = Pm - pgm                # allele residuals
        for i in range(len(tv)):
            if np.std(Pm[i]) > 1e-9 and np.std(T[i]) > 1e-9:
                full_p.append(pearsonr(Pm[i], T[i])[0])
            if np.std(rP[i]) > 1e-9 and np.std(rT[i]) > 1e-9:
                resid_p.append(pearsonr(rP[i], rT[i])[0])
    if full_p:
        rowsB.append(dict(method=m, full_pearson=round(float(np.mean(full_p)), 3),
                          residual_pearson=round(float(np.mean(resid_p)), 3) if resid_p else np.nan,
                          n=len(full_p)))
dfB = pd.DataFrame(rowsB).sort_values('full_pearson', ascending=False)
print(dfB.to_string(index=False))
dfB.to_csv(OUT_DIR / "residual_decomp_models.csv", index=False)
print("\nsaved -> " + str(OUT_DIR) + "/residual_decomp_{measurement,models}.csv")
