#!/usr/bin/env python3
"""Top-3 #3 (structure axis): AlphaMissense as a purpose-built structure+evolution VEP control.

(1) Does AM pathogenicity capture the MEASURED effect magnitude?  Spearman(AM, ||delta_v||).
(2) Can AM rank the transcriptional profiles?  Feed AM (and AM+theta) to the same 6 heads ->
    PDS-cosine + Pearson-delta (expected: PDS ~0.5, since a pathogenicity scalar cannot recover
    the allele-specific expression profile that is floored by measurement).
Matches variant names to AM's protein_variant (WT-pos-MUT); GATA1 multi-substitution variants
have no AM entry and are skipped (reported).
"""
import numpy as np, pandas as pd, sys, os
from scipy.stats import spearmanr, pearsonr
from sklearn.linear_model import Ridge, Lasso
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.decomposition import PCA
sys.path.insert(0, "/data/boom/NUS/VCCompass/unified")
import harness as H

BASE = "/data/boom/NUS/VCCompass"; OUT = f"{BASE}/esm2_control"
UNIPROT = {'TP53': 'P04637', 'KRAS': 'P01116', 'GATA1': 'P15976', 'JAK1': 'P23458'}
U2G = {v: k for k, v in UNIPROT.items()}
SPLITS = ['split1', 'split2', 'split3', 'split5', 'split6']
real = dict(np.load(f"{BASE}/unified/real_deltas.npz"))

am = pd.read_csv(f"{OUT}/AM_4genes.tsv", sep='\t', header=None,
                 names=['uniprot', 'variant', 'am_path', 'am_class'])
amscore = {}
for _, r in am.iterrows():
    g = U2G.get(str(r['uniprot']).strip())
    if g:
        amscore[f"{g}__{str(r['variant']).strip()}"] = float(r['am_path'])
print(f"AlphaMissense scores matched: {len(amscore)}")

# (1) AM pathogenicity vs measured effect size
print("\n=== (1) Spearman(AM pathogenicity, measured ||delta_v||) ===")
allx, ally = [], []
for g in H.GENES:
    xs = [amscore[k] for k in real if k.startswith(g + "__") and k in amscore]
    ys = [np.linalg.norm(real[k]) for k in real if k.startswith(g + "__") and k in amscore]
    if len(xs) >= 5:
        print(f"  {g}: n={len(xs)}  Spearman={spearmanr(xs, ys)[0]:.3f}")
        allx += xs; ally += ys
print(f"  POOLED: n={len(allx)}  Spearman={spearmanr(allx, ally)[0]:.3f}")


def head(name, Xtr, Ytr, Xte):
    G = Ytr.shape[1]
    if name in ('RF', 'GBoost') and G > 200:
        q = min(50, Ytr.shape[0] - 1, G); pca = PCA(q, random_state=0).fit(Ytr)
        base = RandomForestRegressor(50, max_depth=6, random_state=0) if name == 'RF' \
            else MultiOutputRegressor(GradientBoostingRegressor(n_estimators=50, max_depth=3, random_state=0))
        base.fit(Xtr, pca.transform(Ytr)); return pca.inverse_transform(base.predict(Xte))
    m = {'Ridge': Ridge(1.0), 'Lasso': Lasso(0.01, max_iter=2000),
         'RF': RandomForestRegressor(50, max_depth=6, random_state=0),
         'GBoost': MultiOutputRegressor(GradientBoostingRegressor(n_estimators=50, max_depth=3, random_state=0)),
         'KNN': KNeighborsRegressor(3), 'MLP': MLPRegressor((128, 64), max_iter=500, random_state=0)}[name]
    m.fit(Xtr, Ytr); return m.predict(Xte)


def pds(pred, ci, cand):
    if np.linalg.norm(pred) < 1e-12: return 0.5
    d = np.array([1 - np.dot(pred, r) / (np.linalg.norm(pred) * np.linalg.norm(r) + 1e-12) for r in cand])
    td = d[ci]; less = int((d < td - 1e-12).sum()); eq = int((np.abs(d - td) <= 1e-12).sum()); n = len(d)
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1) if n > 1 else 0.5


HEADS = ['Ridge', 'Lasso', 'RF', 'GBoost', 'KNN', 'MLP']
theta = {f"{g}__{v}": np.asarray(vec, np.float32) for g in H.GENES for v, vec in H.theta_map(g).items()}
print("\n=== (2) AM as a feature through the 6 heads: PDS-cosine + Pearson-delta ===")
for featname in ['alphamissense', 'alphamissense+theta']:
    if featname == 'alphamissense':
        feat = {k: np.array([amscore[k]], np.float32) for k in amscore}
    else:
        feat = {k: np.concatenate([[amscore[k]], theta[k]]).astype(np.float32) for k in amscore if k in theta}
    P, E = [], []
    for g in H.GENES:
        for s in SPLITS:
            tr, te = H.split_vars(g, s)
            tr = [v for v in tr if f"{g}__{v}" in real and f"{g}__{v}" in feat]
            cand = [v for v in (tr + te) if f"{g}__{v}" in real and f"{g}__{v}" in feat]
            te = [v for v in te if v in cand]
            if len(tr) < 3 or len(te) < 1: continue
            Xtr = np.stack([feat[f"{g}__{v}"] for v in tr]); Ytr = np.stack([real[f"{g}__{v}"] for v in tr])
            Xte = np.stack([feat[f"{g}__{v}"] for v in te]); Cd = [real[f"{g}__{v}"] for v in cand]
            cidx = {v: i for i, v in enumerate(cand)}
            for h in HEADS:
                try:
                    Yh = head(h, Xtr, Ytr, Xte)
                except Exception:
                    continue
                for i, v in enumerate(te):
                    P.append(pds(Yh[i], cidx[v], Cd))
                    tv = real[f"{g}__{v}"]
                    if np.std(Yh[i]) > 1e-9 and np.std(tv) > 1e-9:
                        E.append(pearsonr(Yh[i], tv)[0])
    print(f"  {featname:22s} dim={len(next(iter(feat.values())))}  PDS={np.mean(P):.3f}  Pearson-delta={np.mean(E):.3f}")
