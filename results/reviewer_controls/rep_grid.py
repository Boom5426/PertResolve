#!/usr/bin/env python3
"""Top-3 #3 scorer: run the SAME 6 regression heads / 5 splits on each representation and
report per-representation PDS-cosine and Pearson-delta (mean over heads). Isolates the
representation variable: identical splits, heads, deterministic pseudobulk, PDS and Pearson.

Revision history that matters for reproducibility:
  The MLP head was constructed as ``MLPRegressor((128, 64), ...)``. In scikit-learn >= 1.7 the
  first positional parameter of MLPRegressor is ``loss``, not ``hidden_layer_sizes``, so every
  MLP fit raised InvalidParameterError and was discarded by a bare ``except Exception``. The
  representation grid therefore ran on FIVE heads, not six. Both defects are fixed here: the
  parameter is passed by keyword, and head-fit failures are counted and reported instead of
  being swallowed. Re-run the grid after this change before citing per-representation numbers.

Usage:
  python rep_grid.py [--base /path/to/VCCompass] [--out /path/to/output_dir]
  VCCOMPASS_BASE=/path/to/VCCompass python rep_grid.py
"""
import argparse
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import Lasso, Ridge
from sklearn.multioutput import MultiOutputRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPRegressor

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_BASE = "/data/boom/NUS/VCCompass"

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--base", default=os.environ.get("VCCOMPASS_BASE", DEFAULT_BASE),
                help="VCCompass compute workspace holding unified/ and esm2_control/ "
                     "(env: VCCOMPASS_BASE)")
ap.add_argument("--out", default=HERE, help="directory to write rep_grid_summary.csv into")
args = ap.parse_args()

BASE = os.path.abspath(args.base)
OUT_DATA = os.path.join(BASE, "esm2_control")
OUT_DIR = os.path.abspath(args.out)
for p in (os.path.join(BASE, "unified", "harness.py"),
          os.path.join(BASE, "unified", "real_deltas.npz"),
          os.path.join(BASE, "esm1v_embeddings.npz")):
    if not os.path.exists(p):
        sys.exit(f"missing required input: {p}\n"
                 f"pass --base or set VCCOMPASS_BASE to the VCCompass workspace")
os.makedirs(OUT_DIR, exist_ok=True)

sys.path.insert(0, os.path.join(BASE, "unified"))
import harness as H  # noqa: E402

SPLITS = ['split1', 'split2', 'split3', 'split5', 'split6']
real = dict(np.load(os.path.join(BASE, "unified", "real_deltas.npz")))
esm1v = dict(np.load(os.path.join(BASE, "esm1v_embeddings.npz")))


def load_rep(name):
    if name == 'theta':
        d = {}
        for g in H.GENES:
            for v, vec in H.theta_map(g).items():
                d[f"{g}__{v}"] = np.asarray(vec, np.float32)
        return d
    if name == 'esm1v':
        return esm1v
    return dict(np.load(os.path.join(OUT_DATA, f"{name}.npz")))


def head(name, Xtr, Ytr, Xte):
    G = Ytr.shape[1]
    if name in ('RF', 'GBoost') and G > 200:
        q = min(50, Ytr.shape[0] - 1, G)
        pca = PCA(q, random_state=0).fit(Ytr)
        Yp = pca.transform(Ytr)
        base = RandomForestRegressor(50, max_depth=6, random_state=0) if name == 'RF' \
            else MultiOutputRegressor(GradientBoostingRegressor(n_estimators=50, max_depth=3, random_state=0))
        base.fit(Xtr, Yp)
        return pca.inverse_transform(base.predict(Xte))
    if name == 'Ridge':
        m = Ridge(alpha=1.0)
    elif name == 'Lasso':
        m = Lasso(alpha=0.01, max_iter=2000)
    elif name == 'RF':
        m = RandomForestRegressor(50, max_depth=6, random_state=0)
    elif name == 'GBoost':
        m = MultiOutputRegressor(GradientBoostingRegressor(n_estimators=50, max_depth=3, random_state=0))
    elif name == 'KNN':
        m = KNeighborsRegressor(n_neighbors=3)
    else:
        # keyword is mandatory: see the module docstring
        m = MLPRegressor(hidden_layer_sizes=(128, 64), max_iter=500, random_state=0)
    m.fit(Xtr, Ytr)
    return m.predict(Xte)


def pds(pred, ci, cand_deltas):
    if np.linalg.norm(pred) < 1e-12:
        return 0.5
    d = np.array([1 - np.dot(pred, r) / (np.linalg.norm(pred) * np.linalg.norm(r) + 1e-12)
                  for r in cand_deltas])
    td = d[ci]
    less = int((d < td - 1e-12).sum())
    eq = int((np.abs(d - td) <= 1e-12).sum())
    n = len(d)
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1) if n > 1 else 0.5


HEADS = ['Ridge', 'Lasso', 'RF', 'GBoost', 'KNN', 'MLP']
REPS = ['theta', 'esm1v', 'esm2_global', 'esm2_globaldelta', 'esm2_sitedelta',
        'esm2_window16', 'esm2_sitedelta+theta']
rows = []
theta = load_rep('theta')
for rep in REPS:
    if rep == 'esm2_sitedelta+theta':
        base_rep = load_rep('esm2_sitedelta')
        feat = {k: np.concatenate([base_rep[k], theta[k]]) for k in base_rep if k in theta}
    else:
        feat = load_rep(rep)
    pds_h, pear_h = {h: [] for h in HEADS}, {h: [] for h in HEADS}
    failures = Counter()
    skipped_cells = 0
    for g in H.GENES:
        for s in SPLITS:
            tr, te = H.split_vars(g, s)
            tr = [v for v in tr if f"{g}__{v}" in real and f"{g}__{v}" in feat]
            cand = [v for v in (tr + te) if f"{g}__{v}" in real and f"{g}__{v}" in feat]
            te = [v for v in te if v in cand]
            if len(tr) < 3 or len(te) < 1:
                skipped_cells += 1
                continue
            Xtr = np.stack([feat[f"{g}__{v}"] for v in tr])
            Ytr = np.stack([real[f"{g}__{v}"] for v in tr])
            Xte = np.stack([feat[f"{g}__{v}"] for v in te])
            Cd = [real[f"{g}__{v}"] for v in cand]
            cidx = {v: i for i, v in enumerate(cand)}
            for h in HEADS:
                try:
                    Yh = head(h, Xtr, Ytr, Xte)
                except Exception as exc:
                    failures[f"{g}/{s}/{h}: {type(exc).__name__}"] += 1
                    continue
                for i, v in enumerate(te):
                    pds_h[h].append(pds(Yh[i], cidx[v], Cd))
                    tv = real[f"{g}__{v}"]
                    if np.std(Yh[i]) > 1e-9 and np.std(tv) > 1e-9:
                        pear_h[h].append(pearsonr(Yh[i], tv)[0])
    heads_used = sorted(h for h in HEADS if pds_h[h])
    allp = [x for h in HEADS for x in pds_h[h]]
    alle = [x for h in HEADS for x in pear_h[h]]
    rows.append(dict(representation=rep, dim=len(next(iter(feat.values()))),
                     PDS_cos=round(float(np.mean(allp)), 3),
                     pearson_delta=round(float(np.mean(alle)), 3),
                     n_heads_scored=len(heads_used), heads=';'.join(heads_used),
                     n_pds=len(allp), n_pearson=len(alle),
                     n_head_failures=int(sum(failures.values()))))
    print(rows[-1], flush=True)
    if failures:
        for k, n in sorted(failures.items()):
            print(f"    head-fit failure {k} x{n}", flush=True)
    if len(heads_used) != len(HEADS):
        print(f"    WARNING: {rep} scored on {len(heads_used)} of {len(HEADS)} heads "
              f"({','.join(heads_used)})", flush=True)

df = pd.DataFrame(rows)
dest = os.path.join(OUT_DIR, "rep_grid_summary.csv")
df.to_csv(dest, index=False)
print("\n" + df.to_string(index=False))
print(f"\nsaved -> {dest}")
