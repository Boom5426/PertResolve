#!/usr/bin/env python3
"""Top-3 #3 scorer: run the SAME 6 regression heads / 5 splits on each representation and
report per-representation PDS-cosine and Pearson-delta (mean over heads). Isolates the
representation variable: identical splits, heads, deterministic pseudobulk, PDS and Pearson.
"""
import numpy as np, pandas as pd, json, sys, os
from scipy.stats import pearsonr
from sklearn.linear_model import Ridge, Lasso
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.decomposition import PCA
sys.path.insert(0, "/data/boom/NUS/VCCompass/unified")
import harness as H

BASE = "/data/boom/NUS/VCCompass"; OUT = f"{BASE}/esm2_control"
SPLITS = ['split1', 'split2', 'split3', 'split5', 'split6']
real = dict(np.load(f"{BASE}/unified/real_deltas.npz"))
esm1v = dict(np.load(f"{BASE}/esm1v_embeddings.npz"))


def load_rep(name):
    if name == 'theta':
        d = {}
        for g in H.GENES:
            for v, vec in H.theta_map(g).items():
                d[f"{g}__{v}"] = np.asarray(vec, np.float32)
        return d
    if name == 'esm1v':
        return esm1v
    return dict(np.load(f"{OUT}/{name}.npz"))


def head(name, Xtr, Ytr, Xte):
    G = Ytr.shape[1]
    if name in ('RF', 'GBoost') and G > 200:
        q = min(50, Ytr.shape[0] - 1, G)
        pca = PCA(q, random_state=0).fit(Ytr); Yp = pca.transform(Ytr)
        base = RandomForestRegressor(50, max_depth=6, random_state=0) if name == 'RF' \
            else MultiOutputRegressor(GradientBoostingRegressor(n_estimators=50, max_depth=3, random_state=0))
        base.fit(Xtr, Yp)
        return pca.inverse_transform(base.predict(Xte))
    if name == 'Ridge': m = Ridge(alpha=1.0)
    elif name == 'Lasso': m = Lasso(alpha=0.01, max_iter=2000)
    elif name == 'RF': m = RandomForestRegressor(50, max_depth=6, random_state=0)
    elif name == 'GBoost': m = MultiOutputRegressor(GradientBoostingRegressor(n_estimators=50, max_depth=3, random_state=0))
    elif name == 'KNN': m = KNeighborsRegressor(n_neighbors=3)
    else: m = MLPRegressor((128, 64), max_iter=500, random_state=0)
    m.fit(Xtr, Ytr)
    return m.predict(Xte)


def pds(pred, ci, cand_deltas):
    if np.linalg.norm(pred) < 1e-12: return 0.5
    d = np.array([1 - np.dot(pred, r) / (np.linalg.norm(pred) * np.linalg.norm(r) + 1e-12) for r in cand_deltas])
    td = d[ci]; less = int((d < td - 1e-12).sum()); eq = int((np.abs(d - td) <= 1e-12).sum())
    n = len(d)
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1) if n > 1 else 0.5


HEADS = ['Ridge', 'Lasso', 'RF', 'GBoost', 'KNN', 'MLP']
REPS = ['theta', 'esm1v', 'esm2_global', 'esm2_globaldelta', 'esm2_sitedelta', 'esm2_window16', 'esm2_sitedelta+theta']
rows = []
theta = load_rep('theta')
for rep in REPS:
    if rep == 'esm2_sitedelta+theta':
        base = load_rep('esm2_sitedelta')
        feat = {k: np.concatenate([base[k], theta[k]]) for k in base if k in theta}
    else:
        feat = load_rep(rep)
    pds_h, pear_h = {h: [] for h in HEADS}, {h: [] for h in HEADS}
    for g in H.GENES:
        for s in SPLITS:
            tr, te = H.split_vars(g, s)
            tr = [v for v in tr if f"{g}__{v}" in real and f"{g}__{v}" in feat]
            cand = [v for v in (tr + te) if f"{g}__{v}" in real and f"{g}__{v}" in feat]
            te = [v for v in te if v in cand]
            if len(tr) < 3 or len(te) < 1: continue
            Xtr = np.stack([feat[f"{g}__{v}"] for v in tr]); Ytr = np.stack([real[f"{g}__{v}"] for v in tr])
            Xte = np.stack([feat[f"{g}__{v}"] for v in te])
            Cd = [real[f"{g}__{v}"] for v in cand]; cidx = {v: i for i, v in enumerate(cand)}
            for h in HEADS:
                try:
                    Yh = head(h, Xtr, Ytr, Xte)
                except Exception:
                    continue
                for i, v in enumerate(te):
                    pds_h[h].append(pds(Yh[i], cidx[v], Cd))
                    tv = real[f"{g}__{v}"]
                    if np.std(Yh[i]) > 1e-9 and np.std(tv) > 1e-9:
                        pear_h[h].append(pearsonr(Yh[i], tv)[0])
    allp = [x for h in HEADS for x in pds_h[h]]; alle = [x for h in HEADS for x in pear_h[h]]
    rows.append(dict(representation=rep, dim=len(next(iter(feat.values()))),
                     PDS_cos=round(float(np.mean(allp)), 3), pearson_delta=round(float(np.mean(alle)), 3)))
    print(rows[-1], flush=True)

df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/rep_grid_summary.csv", index=False)
print("\n" + df.to_string(index=False))
print(f"\nsaved -> {OUT}/rep_grid_summary.csv")
