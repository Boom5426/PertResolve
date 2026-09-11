#!/usr/bin/env python3
"""Decisive test: does confining predictions to PerturbNet's 50-dim WT-PCA subspace
inflate PDS? Project a KNOWN-at-chance predictor (Ridge-esm, full-space PDS ~0.476)
into that subspace and re-score. Also a pure random-in-subspace baseline. If these
jump to ~0.56, PerturbNet's 0.56 is a subspace artifact, not prediction.

The six PDS values are printed exactly as before and additionally written, at full
precision, to the JSON file named by --out. Those are the six numbers quoted in
results/canonical/perturbnet_subspace_test.json; the permutation-null fields in that
canonical table come from a separate script and are not recomputed here.

Usage:
  python subspace_test.py --out /path/to/subspace_test.json [--base /path/to/processed-data]
  PERTRESOLVE_DATA=/path/to/processed-data python subspace_test.py --out /path/to/subspace_test.json
"""
import numpy as np, sys, glob, os
import argparse, json
from pathlib import Path

# The package is not installed by default; make this repository importable so that
# the shared path helper can be used when the script is run directly.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import (resolve_base, add_harness_to_path, require_inputs,
                                 reject_repo_results)

_ap = argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--base", default=None,
                 help="directory holding the shared scorer (env: PERTRESOLVE_DATA)")
_ap.add_argument("--out", required=True,
                 help="path of the JSON file to write the six PDS values into")
_args = _ap.parse_args()

BASE = resolve_base(_args.base)
OUT_PATH = reject_repo_results(_args.out)
require_inputs(BASE / "unified" / "real_deltas.npz",
               BASE / "unified" / "preds5" / "Ridge-esm.npz",
               BASE / "unified" / "preds5" / "PerturbNet.npz",
               BASE / "unified" / "preds5" / "Gene-mean.npz")
add_harness_to_path(BASE)

import harness as H
from sklearn.decomposition import PCA
np.random.seed(0)
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

ridge_esm_full_pds = score(re_full)
print(f"Ridge-esm  FULL space PDS      = {ridge_esm_full_pds:.3f}   (reference: ~0.476, at chance)")
ridge_esm_subspace_pds = score(re_sub)
print(f"Ridge-esm  PROJECTED to subspace = {ridge_esm_subspace_pds:.3f}")
gene_mean_full_pds = score(gm_full)
print(f"Gene-mean  FULL space PDS      = {gene_mean_full_pds:.3f}")
gene_mean_subspace_pds = score(gm_sub)
print(f"Gene-mean  PROJECTED to subspace = {gene_mean_subspace_pds:.3f}")
pure_random_in_subspace_pds = score(rand_sub)
print(f"PURE RANDOM in subspace PDS    = {pure_random_in_subspace_pds:.3f}   (any value >0.5 => pure artifact)")
perturbnet_subspace_pds = score(pn)
print(f"PerturbNet (as scored)         = {perturbnet_subspace_pds:.3f}")

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
OUT_PATH.write_text(json.dumps({
    "ridge_esm_full_pds": ridge_esm_full_pds,
    "ridge_esm_subspace_pds": ridge_esm_subspace_pds,
    "gene_mean_full_pds": gene_mean_full_pds,
    "gene_mean_subspace_pds": gene_mean_subspace_pds,
    "pure_random_in_subspace_pds": pure_random_in_subspace_pds,
    "perturbnet_subspace_pds": perturbnet_subspace_pds,
}, indent=2) + "\n")
