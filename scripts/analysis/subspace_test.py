#!/usr/bin/env python3
"""Decisive test: does confining predictions to PerturbNet's 50-dim WT-PCA subspace
inflate PDS? Project a KNOWN-at-chance predictor (Ridge-esm, full-space PDS ~0.476)
into that subspace and re-score. Also a pure random-in-subspace baseline. If these
jump to ~0.56, PerturbNet's 0.56 is a subspace artifact, not prediction."""
import numpy as np, sys, glob, os
sys.path.insert(0, "/data/boom/NUS/VCCompass/unified"); import harness as H
from sklearn.decomposition import PCA
np.random.seed(0); BASE = "/data/boom/NUS/VCCompass"
SPLITS = ['split1','split2','split3','split5','split6']
real = {k: v for k, v in np.load(f"{BASE}/unified/real_deltas.npz").items()}
cand = {(g, s): [v for v in (H.split_vars(g, s)[0] + H.split_vars(g, s)[1]) if f"{g}__{v}" in real] for g in H.GENES for s in SPLITS}

# per-gene WT-PCA-50 components (identical to perturbnet_canonical.py)
comps = {}
for g in H.GENES:
    X, lab = H.load_gene(g); lab = np.asarray(lab); wm = lab == 'WT'
    npc = min(50, X.shape[1] - 1)
    comps[g] = PCA(npc, random_state=42).fit(X[wm][:min(5000, wm.sum())]).components_.astype(np.float32)

def proj(g, x):
    C = comps[g]; return (x @ C.T) @ C   # orthogonal projection onto the 50-dim subspace

def score(pred_lookup):  # pred_lookup(g,s,v) -> delta ; returns overall mean-of-(split,gene)
    cells = []
    for g in H.GENES:
        for s in SPLITS:
            te = H.split_vars(g, s)[1]; te = [v for v in te if f"{g}__{v}" in real and v in cand[(g,s)]]
            cv = cand[(g, s)]; rd = {u: real[f"{g}__{u}"] for u in cv}
            pv = []
            for v in te:
                p = pred_lookup(g, s, v)
                if p is None: continue
                pv.append(H.canonical_pds_cos(p, v, cv, rd))
            if pv: cells.append(np.mean(pv))
    return float(np.mean(cells))

# Ridge-esm predictions (full space, known ~0.476)
RE = np.load(f"{BASE}/unified/preds5/Ridge-esm.npz")
def re_full(g, s, v): k = f"{g}__{s}__{v}"; return RE[k] if k in RE.files else None
def re_sub(g, s, v):  k = f"{g}__{s}__{v}"; return proj(g, RE[k]) if k in RE.files else None
# PerturbNet (already subspace-reconstructed)
PN = np.load(f"{BASE}/unified/preds5/PerturbNet.npz")
def pn(g, s, v): k = f"{g}__{s}__{v}"; return PN[k] if k in PN.files else None
# pure random-in-subspace: random 50-dim -> gene space
rng = np.random.RandomState(1)
randcache = {}
def rand_sub(g, s, v):
    key = (g, s, v)
    if key not in randcache:
        C = comps[g]; z = rng.randn(C.shape[0]).astype(np.float32); randcache[key] = z @ C
    return randcache[key]
# Gene-mean full vs sub
GM = np.load(f"{BASE}/unified/preds5/Gene-mean.npz")
def gm_full(g, s, v): k = f"{g}__{s}__{v}"; return GM[k] if k in GM.files else None
def gm_sub(g, s, v):  k = f"{g}__{s}__{v}"; return proj(g, GM[k]) if k in GM.files else None

print(f"Ridge-esm  FULL space PDS      = {score(re_full):.3f}   (reference: ~0.476, at chance)")
print(f"Ridge-esm  PROJECTED to subspace = {score(re_sub):.3f}")
print(f"Gene-mean  FULL space PDS      = {score(gm_full):.3f}")
print(f"Gene-mean  PROJECTED to subspace = {score(gm_sub):.3f}")
print(f"PURE RANDOM in subspace PDS    = {score(rand_sub):.3f}   (any value >0.5 => pure artifact)")
print(f"PerturbNet (as scored)         = {score(pn):.3f}")
