#!/usr/bin/env python3
"""Top-3 #3 (structure axis): AlphaMissense as a purpose-built structure+evolution VEP control.

(1) Does AM pathogenicity capture the MEASURED effect magnitude?  Spearman(AM, ||delta_v||).
(2) Can AM rank the transcriptional profiles?  Feed AM (and AM+theta) to the same 6 heads ->
    PDS-cosine + Pearson-delta (expected: PDS ~0.5, since a pathogenicity scalar cannot recover
    the allele-specific expression profile that is floored by measurement).
Matches variant names to AM's protein_variant (WT-pos-MUT); GATA1 multi-substitution variants
have no AM entry and are skipped (reported).

Outputs (written next to this script unless --out is given):
  am_summary.csv          representation, dim, PDS_cos, pearson_delta, n_pds, n_pearson,
                          n_head_failures  -- schema-compatible with rep_grid_summary.csv
  am_effect_spearman.csv  gene, n, spearman

The scoring logic is unchanged from the run that produced the manuscript values; this revision
only removes hard-coded paths, writes the results to tracked files instead of printing them,
and reports head-fit failures explicitly instead of swallowing them.

Usage:
  python am_analysis.py [--base /path/to/processed-data] [--out /path/to/output_dir]
  PERTRESOLVE_DATA=/path/to/processed-data python am_analysis.py
"""
import argparse
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import Lasso, Ridge
from sklearn.multioutput import MultiOutputRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPRegressor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
from pertresolve.paths import reject_repo_results, resolve_base  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--base", default=None,
                help="processed-data workspace holding unified/ and esm2_control/ "
                     "(env: PERTRESOLVE_DATA)")
# This used to default to HERE, which is inside results/: a bare re-run overwrote the
# committed tables in place, which is exactly what reject_repo_results exists to stop.
# Naming a scratch directory is now required, and promotion stays a separate, deliberate step.
ap.add_argument("--out", required=True,
                help="directory to write the result tables into; may not be inside results/")
args = ap.parse_args()

BASE = str(resolve_base(args.base, what="the processed-data workspace"))
OUT_DATA = os.path.join(BASE, "esm2_control")
OUT_DIR = str(reject_repo_results(args.out))
for p in (os.path.join(BASE, "unified", "harness.py"),
          os.path.join(BASE, "unified", "real_deltas.npz"),
          os.path.join(OUT_DATA, "AM_4genes.tsv")):
    if not os.path.exists(p):
        sys.exit(f"missing required input: {p}\n"
                 f"pass --base or set PERTRESOLVE_DATA to the processed-data workspace")
os.makedirs(OUT_DIR, exist_ok=True)

sys.path.insert(0, os.path.join(BASE, "unified"))
import harness as H  # noqa: E402

UNIPROT = {'TP53': 'P04637', 'KRAS': 'P01116', 'GATA1': 'P15976', 'JAK1': 'P23458'}
U2G = {v: k for k, v in UNIPROT.items()}
SPLITS = ['split1', 'split2', 'split3', 'split5', 'split6']
real = dict(np.load(os.path.join(BASE, "unified", "real_deltas.npz")))

am = pd.read_csv(os.path.join(OUT_DATA, "AM_4genes.tsv"), sep='\t', header=None,
                 names=['uniprot', 'variant', 'am_path', 'am_class'])
amscore = {}
for _, r in am.iterrows():
    g = U2G.get(str(r['uniprot']).strip())
    if g:
        amscore[f"{g}__{str(r['variant']).strip()}"] = float(r['am_path'])
print(f"AlphaMissense rows read: {len(am)};  keys matched to the four genes: {len(amscore)}")

# Coverage: which measured variants have no AM entry (GATA1 multi-substitution alleles).
print("\n=== AM coverage of the measured variant set ===")
coverage = []
for g in H.GENES:
    meas = [k for k in real if k.startswith(g + "__")]
    hit = [k for k in meas if k in amscore]
    coverage.append((g, len(meas), len(hit)))
    print(f"  {g}: {len(hit)}/{len(meas)} measured variants have an AM score")

# (1) AM pathogenicity vs measured effect size
print("\n=== (1) Spearman(AM pathogenicity, measured ||delta_v||) ===")
allx, ally, sp_rows = [], [], []
for g in H.GENES:
    xs = [amscore[k] for k in real if k.startswith(g + "__") and k in amscore]
    ys = [np.linalg.norm(real[k]) for k in real if k.startswith(g + "__") and k in amscore]
    if len(xs) >= 5:
        rho = spearmanr(xs, ys)[0]
        print(f"  {g}: n={len(xs)}  Spearman={rho:.3f}")
        sp_rows.append({'gene': g, 'n': len(xs), 'spearman': round(float(rho), 4)})
        allx += xs
        ally += ys
    else:
        print(f"  {g}: n={len(xs)}  (fewer than 5 AM-scored variants, not correlated)")
        sp_rows.append({'gene': g, 'n': len(xs), 'spearman': np.nan})
pooled = spearmanr(allx, ally)[0]
print(f"  POOLED: n={len(allx)}  Spearman={pooled:.3f}")
sp_rows.append({'gene': 'POOLED', 'n': len(allx), 'spearman': round(float(pooled), 4)})


def head(name, Xtr, Ytr, Xte):
    G = Ytr.shape[1]
    if name in ('RF', 'GBoost') and G > 200:
        q = min(50, Ytr.shape[0] - 1, G)
        pca = PCA(q, random_state=0).fit(Ytr)
        base = RandomForestRegressor(50, max_depth=6, random_state=0) if name == 'RF' \
            else MultiOutputRegressor(GradientBoostingRegressor(n_estimators=50, max_depth=3, random_state=0))
        base.fit(Xtr, pca.transform(Ytr))
        return pca.inverse_transform(base.predict(Xte))
    m = {'Ridge': Ridge(1.0), 'Lasso': Lasso(0.01, max_iter=2000),
         'RF': RandomForestRegressor(50, max_depth=6, random_state=0),
         'GBoost': MultiOutputRegressor(GradientBoostingRegressor(n_estimators=50, max_depth=3, random_state=0)),
         'KNN': KNeighborsRegressor(3),
         # hidden_layer_sizes MUST be keyword: MLPRegressor's first positional parameter is
         # `loss` in scikit-learn >= 1.7, so the old positional form raised
         # InvalidParameterError on every fit and the MLP head silently dropped out.
         'MLP': MLPRegressor(hidden_layer_sizes=(128, 64), max_iter=500, random_state=0)}[name]
    m.fit(Xtr, Ytr)
    return m.predict(Xte)


def pds(pred, ci, cand):
    if np.linalg.norm(pred) < 1e-12:
        return 0.5
    d = np.array([1 - np.dot(pred, r) / (np.linalg.norm(pred) * np.linalg.norm(r) + 1e-12) for r in cand])
    td = d[ci]
    less = int((d < td - 1e-12).sum())
    eq = int((np.abs(d - td) <= 1e-12).sum())
    n = len(d)
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1) if n > 1 else 0.5


HEADS = ['Ridge', 'Lasso', 'RF', 'GBoost', 'KNN', 'MLP']
theta = {f"{g}__{v}": np.asarray(vec, np.float32) for g in H.GENES for v, vec in H.theta_map(g).items()}
print("\n=== (2) AM as a feature through the 6 heads: PDS-cosine + Pearson-delta ===")
grid_rows = []
for featname in ['alphamissense', 'alphamissense+theta']:
    if featname == 'alphamissense':
        feat = {k: np.array([amscore[k]], np.float32) for k in amscore}
    else:
        feat = {k: np.concatenate([[amscore[k]], theta[k]]).astype(np.float32)
                for k in amscore if k in theta}
    P, E = [], []
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
                    # Control flow is deliberately unchanged from the manuscript run; the
                    # failure is now counted and reported rather than silently dropped.
                    failures[f"{g}/{s}/{h}: {type(exc).__name__}"] += 1
                    continue
                for i, v in enumerate(te):
                    P.append(pds(Yh[i], cidx[v], Cd))
                    tv = real[f"{g}__{v}"]
                    if np.std(Yh[i]) > 1e-9 and np.std(tv) > 1e-9:
                        E.append(pearsonr(Yh[i], tv)[0])
    dim = len(next(iter(feat.values())))
    print(f"  {featname:22s} dim={dim}  PDS={np.mean(P):.3f}  Pearson-delta={np.mean(E):.3f}"
          f"  (n_pds={len(P)}, n_pearson={len(E)}, gene-split cells skipped for size={skipped_cells},"
          f" head fits failed={sum(failures.values())})")
    if failures:
        print("     head-fit failures:")
        for k, n in sorted(failures.items()):
            print(f"       {k} x{n}")
    grid_rows.append({'representation': featname, 'dim': dim,
                      'PDS_cos': round(float(np.mean(P)), 3),
                      'pearson_delta': round(float(np.mean(E)), 3),
                      'n_pds': len(P), 'n_pearson': len(E),
                      'n_head_failures': int(sum(failures.values()))})

grid = os.path.join(OUT_DIR, "am_summary.csv")
sp = os.path.join(OUT_DIR, "am_effect_spearman.csv")
pd.DataFrame(grid_rows).to_csv(grid, index=False)
pd.DataFrame(sp_rows).to_csv(sp, index=False)
print(f"\nwrote {grid}")
print(f"wrote {sp}")
print("\nprovenance: AM scores from AlphaMissense_aa_substitutions.tsv.gz filtered to the four "
      "UniProt accessions (P04637 TP53, P01116 KRAS, P15976 GATA1, P23458 JAK1); measured "
      "targets from unified/real_deltas.npz; splits and theta from unified/harness.py. All "
      "estimators use random_state=0.")
