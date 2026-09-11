#!/usr/bin/env python3
"""1D: analytical scaling law g(rho) and collapse.

Derivation (Euclidean pseudobulk means, stated in Methods): with per-profile sampling
noise eta^2 = tr(Sigma)/m, two disjoint half-profiles of the SAME variant satisfy
  E||d_A - d_B||^2 = 2 eta^2,
and two different variants satisfy
  E||d_i - d_j||^2 = ||mu_i - mu_j||^2 + 2 eta^2.
So the between-variant signal Delta^2 = mean_pair ||d_i - d_j||^2 - 2 eta^2 (noise-debiased),
and the dimensionless discrimination SNR is
  rho = sqrt(Delta^2 / (2 eta^2)).
Because eta ~ 1/sqrt(m), rho carries the sqrt(n) depth scaling. Claim: the oracle ceiling
(and PDS) is a single monotone function g(rho), so all (gene, depth) points collapse onto
one curve; g is calibrated here, its scaling is derived. rho and eta are estimated from
disjoint halves (no circularity with the cosine PDS ceiling they are compared against).
"""
import numpy as np, pandas as pd, sys
from scipy.stats import spearmanr
import argparse
from pathlib import Path

# The package is not installed by default; make this repository importable so that
# the shared path helpers can be used when the script is run directly.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import (
    add_harness_to_path,
    reject_repo_results,
    require_inputs,
    resolve_base,
)

_ap = argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--base", default=None,
                 help="directory holding the shared scorer and the gene "
                      "arrays (env: PERTRESOLVE_DATA)")
_ap.add_argument("--out", required=True, metavar="OUT_DIR",
                 help="directory receiving floor_law.csv; may not be inside the "
                      "repository's results/")
_args = _ap.parse_args()

BASE = resolve_base(_args.base)
OUT_DIR = reject_repo_results(_args.out)
OUT_DIR.mkdir(parents=True, exist_ok=True)
add_harness_to_path(BASE)

import harness as H

NSEED = 20
WT_TAGS = ('WT', 'wt', 'WT_control')
DEPTHS = [25, 50, 100, 150, 250, 400]
np.random.seed(0)
gene_cells = {g: (lambda X, lab: (X, np.asarray(lab)))(*H.load_gene(g)) for g in H.GENES}


def pds_row(drow, ci, n):
    td = drow[ci]
    if not np.isfinite(td):
        return 0.5
    less = int((drow < td - 1e-12).sum())
    eq = int((np.abs(drow - td) <= 1e-12).sum())
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1) if n > 1 else 0.5


def norm(M):
    return M / (np.linalg.norm(M, axis=1, keepdims=True) + 1e-12)


def build_eval(g, m, seed):
    X, lab = gene_cells[g]
    rng = np.random.RandomState(1000 + 7 * seed + H._GENE_SEED[g])
    wt = np.where(np.isin(lab, WT_TAGS))[0]
    if len(wt) < 2 * m:
        return None, None
    rng.shuffle(wt)
    wm_b, wm_e = X[wt[:m]].mean(0), X[wt[m:2 * m]].mean(0)
    bd, ed = {}, {}
    for v in np.unique(lab):
        if v in WT_TAGS:
            continue
        idx = np.where(lab == v)[0]
        if len(idx) < 2 * m:
            continue
        rng.shuffle(idx)
        bd[v] = (X[idx[:m]].mean(0) - wm_b).astype(np.float32)
        ed[v] = (X[idx[m:2 * m]].mean(0) - wm_e).astype(np.float32)
    return bd, ed


def mean_pair_sqdist(M):
    # mean over i<j of ||M_i - M_j||^2
    n = M.shape[0]
    G = M @ M.T
    sq = np.diag(G)
    D2 = sq[:, None] + sq[None, :] - 2 * G
    iu = np.triu_indices(n, 1)
    return float(D2[iu].mean())


rows = []
for g in H.GENES:
    for m in DEPTHS:
        eta2_s, pair2_s, ceil_s, n_s = [], [], [], None
        for seed in range(NSEED):
            bd, ed = build_eval(g, m, seed)
            if bd is None or len(bd) < 5:
                continue
            vs = sorted(bd.keys())
            n_s = len(vs)
            B = np.stack([bd[v] for v in vs])
            E = np.stack([ed[v] for v in vs])
            eta2_s.append(0.5 * float(np.mean(np.sum((B - E) ** 2, axis=1))))  # per-profile noise^2
            pair2_s.append(mean_pair_sqdist(E))                                # observed between-variant
            # oracle ceiling (cosine PDS), query=build, truth=eval
            D = 1 - norm(B) @ norm(E).T
            ceil_s.append(float(np.mean([pds_row(D[i], i, n_s) for i in range(n_s)])))
        if n_s is None or not eta2_s:
            continue
        # seed-averaged estimators, THEN debias once (law of large numbers avoids per-seed cancellation)
        eta2 = float(np.mean(eta2_s))
        pair2 = float(np.mean(pair2_s))
        delta2 = max(pair2 - 2 * eta2, 0.0)
        rho = float(np.sqrt(delta2 / (2 * eta2 + 1e-12)))
        rows.append(dict(gene=g, depth_m=m, n_var=n_s,
                         eta2=round(eta2, 2), delta2=round(delta2, 2),
                         rho=round(rho, 3),
                         ceiling_pds=round(float(np.mean(ceil_s)), 3)))

df = pd.DataFrame(rows).sort_values('rho')
df.to_csv(OUT_DIR / "floor_law.csv", index=False)
print(df.to_string(index=False))

sr, _ = spearmanr(df['rho'], df['ceiling_pds'])
print(f"\nSpearman(rho, ceiling) across ALL gene-depth points = {sr:.3f} "
      f"(collapse onto one monotone g(rho) => close to 1)")
# show that different genes at matched rho have matched ceiling (collapse check)
print("\n=== collapse check: bin by rho, spread of ceiling across genes within a bin ===")
df['rho_bin'] = pd.cut(df['rho'], [0, 0.5, 1.0, 1.5, 2.5, 10],
                       labels=['<0.5', '0.5-1', '1-1.5', '1.5-2.5', '>2.5'])
print(df.groupby('rho_bin', observed=True).agg(
    n_points=('ceiling_pds', 'size'),
    genes=('gene', lambda s: ','.join(sorted(set(s)))),
    ceiling_mean=('ceiling_pds', 'mean'),
    ceiling_std=('ceiling_pds', 'std')).round(3).to_string())
print(f"\nsaved -> {OUT_DIR / 'floor_law.csv'}")
