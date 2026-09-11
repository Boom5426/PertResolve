#!/usr/bin/env python3
"""2A part 1: metric-family invariance of the measurement floor.

Recompute the split-half window ratio D_self/D_null under five distance families in the
canonical PCA-50 space, to show the floor is not an artifact of energy distance.
  D_self = dist(halfA, halfB)   [two disjoint halves of the SAME variant]
  D_null = dist(halfA, WT sample)
Distributional metrics (energy, MMD-RBF, sliced-Wasserstein) act on the two cell SETS;
mean metrics (cosine, L2) act on the set means. Per gene we report median(D_self)/median(D_null)
(the manuscript convention). Energy should reproduce the canonical Fig 3b ratios as a check.
"""
import numpy as np, pandas as pd, sys
from sklearn.decomposition import PCA
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
                 help="directory receiving metric_family_floor.csv; may not be inside the "
                      "repository's results/")
_args = _ap.parse_args()

BASE = resolve_base(_args.base)
OUT_DIR = reject_repo_results(_args.out)
OUT_DIR.mkdir(parents=True, exist_ok=True)
add_harness_to_path(BASE)

import harness as H

M = 50            # cells per half
NSEED = 10
MAXVAR = 40       # cap variants per gene for speed (random subset with enough cells)
LPROJ = 100       # slices for sliced-Wasserstein
WT_TAGS = ('WT', 'wt', 'WT_control')
np.random.seed(0)


def energy(A, B):
    dab = np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1)).mean()
    daa = np.sqrt(((A[:, None, :] - A[None, :, :]) ** 2).sum(-1)).mean()
    dbb = np.sqrt(((B[:, None, :] - B[None, :, :]) ** 2).sum(-1)).mean()
    return 2 * dab - daa - dbb


def mmd_rbf(A, B, sig):
    def k(X, Y):
        d2 = ((X[:, None, :] - Y[None, :, :]) ** 2).sum(-1)
        return np.exp(-d2 / (2 * sig * sig))
    return max(k(A, A).mean() + k(B, B).mean() - 2 * k(A, B).mean(), 0.0) ** 0.5


def sliced_w(A, B, rng):
    P = rng.standard_normal((A.shape[1], LPROJ)); P /= np.linalg.norm(P, axis=0, keepdims=True) + 1e-12
    pa = np.sort(A @ P, axis=0); pb = np.sort(B @ P, axis=0)
    return float(np.abs(pa - pb).mean())


def cos_mean(A, B):
    a, b = A.mean(0), B.mean(0)
    return 1 - float(a @ b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12)


def l2_mean(A, B):
    return float(np.linalg.norm(A.mean(0) - B.mean(0)))


rows = []
for g in H.GENES:
    X, lab = H.load_gene(g); lab = np.asarray(lab)
    Xp = PCA(n_components=50, random_state=0).fit_transform(X.astype(np.float32))
    wt_idx = np.where(np.isin(lab, WT_TAGS))[0]
    variants = [v for v in np.unique(lab) if v not in WT_TAGS and (lab == v).sum() >= 2 * M]
    rng0 = np.random.RandomState(0)
    if len(variants) > MAXVAR:
        variants = list(rng0.choice(variants, MAXVAR, replace=False))
    # global sigma (median heuristic) for MMD in this gene's PCA space
    samp = Xp[rng0.choice(len(Xp), min(500, len(Xp)), replace=False)]
    d2 = ((samp[:, None, :] - samp[None, :, :]) ** 2).sum(-1)
    sig = np.sqrt(np.median(d2[d2 > 0]) / 2) + 1e-9
    per = {m: {'self': [], 'null': []} for m in ['energy', 'mmd', 'slicedW', 'cosine', 'L2']}
    for seed in range(NSEED):
        rng = np.random.RandomState(100 + seed + H._GENE_SEED[g])
        wt = wt_idx.copy(); rng.shuffle(wt); C = Xp[wt[:M]] if len(wt) >= M else Xp[wt]
        for v in variants:
            idx = np.where(lab == v)[0]; rng.shuffle(idx)
            A, B = Xp[idx[:M]], Xp[idx[M:2 * M]]
            per['energy']['self'].append(energy(A, B));   per['energy']['null'].append(energy(A, C))
            per['mmd']['self'].append(mmd_rbf(A, B, sig)); per['mmd']['null'].append(mmd_rbf(A, C, sig))
            per['slicedW']['self'].append(sliced_w(A, B, rng)); per['slicedW']['null'].append(sliced_w(A, C, rng))
            per['cosine']['self'].append(cos_mean(A, B)); per['cosine']['null'].append(cos_mean(A, C))
            per['L2']['self'].append(l2_mean(A, B));       per['L2']['null'].append(l2_mean(A, C))
    for m in per:
        rs = np.median(per[m]['self']) / (np.median(per[m]['null']) + 1e-12)
        rows.append(dict(gene=g, metric=m, ratio=round(float(rs), 3),
                         med_Dself=round(float(np.median(per[m]['self'])), 3),
                         med_Dnull=round(float(np.median(per[m]['null'])), 3)))

df = pd.DataFrame(rows)
piv = df.pivot(index='gene', columns='metric', values='ratio').reindex(H.GENES)
piv = piv[['energy', 'mmd', 'slicedW', 'cosine', 'L2']]
df.to_csv(OUT_DIR / "metric_family_floor.csv", index=False)
print("D_self/D_null ratio by gene x metric (low = wide window / rankable; ~1 = at floor):")
print(piv.round(3).to_string())
print(f"\nsaved -> {OUT_DIR / 'metric_family_floor.csv'}")
