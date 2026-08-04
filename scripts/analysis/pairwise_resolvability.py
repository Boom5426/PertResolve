#!/usr/bin/env python3
"""1B: Pairwise allele-allele identification.

PDS measures v-vs-other-variants (identification), not v-vs-WT (detection). A variant
easy to detect vs WT can still be unresolvable from its sibling alleles. All quantities
are computed at a single matched split-half depth in the PDS cosine-delta space:

  d_A(v), d_B(v)  : disjoint split-half pseudobulk deltas (same construction as 1A)
  D_self(v)       : cos_dist(d_A, d_B), seed-averaged   [replicate noise]
  D(v_i, v_j)     : cos_dist(d_A(v_i), d_A(v_j)), seed-averaged   [between-variant signal]

A variant is 'identifiable' if its nearest sibling is beyond its own replicate noise
(min_j D(i,j) > D_self(i)); a pair is 'resolvable' if D(i,j) > mean(D_self_i, D_self_j).
Candidate siblings = all variants of the same gene (identification is a gene property).
"""
import numpy as np, pandas as pd, sys
import argparse
from pathlib import Path

# The package is not installed by default; make this repository importable so that
# the shared path helpers can be used when the script is run directly.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from alleleperturb.paths import (
    add_harness_to_path,
    reject_repo_results,
    require_inputs,
    resolve_base,
)

_ap = argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--base", default=None,
                 help="directory holding the shared scorer and the gene "
                      "arrays (env: ALLELEPERTURB_DATA)")
_ap.add_argument("--out", required=True, metavar="OUT_DIR",
                 help="directory receiving pairwise_resolvability.csv; may not be inside the "
                      "repository's results/")
_args = _ap.parse_args()

BASE = resolve_base(_args.base)
OUT_DIR = reject_repo_results(_args.out)
OUT_DIR.mkdir(parents=True, exist_ok=True)
require_inputs(BASE / "allele_perturb_bench.csv")
add_harness_to_path(BASE)

import harness as H

NSUB = 300
NSEED = 15
WT_TAGS = ('WT', 'wt', 'WT_control')
np.random.seed(0)
gene_cells = {g: (lambda X, lab: (X, np.asarray(lab)))(*H.load_gene(g)) for g in H.GENES}
_BENCH_DF = pd.read_csv(f"{BASE}/allele_perturb_bench.csv")
BENCH = {g: set(_BENCH_DF[_BENCH_DF.gene == g]["variant"]) for g in H.GENES}


def halves(seed, g):
    X, lab = gene_cells[g]
    rng = np.random.RandomState(1000 + 7 * seed + H._GENE_SEED[g])
    wt = np.where(np.isin(lab, WT_TAGS))[0]
    if len(wt) == 0:
        wm_t = wm_p = np.zeros(X.shape[1], np.float32)
    else:
        rng.shuffle(wt)
        h = len(wt) // 2
        a, b = wt[:h][:NSUB], wt[h:][:NSUB]
        wm_t, wm_p = X[a].mean(0), X[b].mean(0)
    dA, dB = {}, {}
    for v in np.unique(lab):
        if v in WT_TAGS or v not in BENCH[g]:
            continue
        idx = np.where(lab == v)[0]
        if len(idx) < 5:
            continue
        rng.shuffle(idx)
        h = len(idx) // 2
        ia, ib = idx[:h][:NSUB], idx[h:][:NSUB]
        dA[v] = (X[ia].mean(0) - wm_t).astype(np.float32)
        dB[v] = (X[ib].mean(0) - wm_p).astype(np.float32)
    return dA, dB


def norm(M):
    return M / (np.linalg.norm(M, axis=1, keepdims=True) + 1e-12)


rows = []
for g in H.GENES:
    variants = None
    dself_acc = {}
    Dpair_acc = None
    for seed in range(NSEED):
        dA, dB = halves(seed, g)
        vs = sorted(dA.keys())
        if variants is None:
            variants = vs
        A = norm(np.stack([dA[v] for v in vs]))
        Dp = 1 - A @ A.T
        Dpair_acc = Dp if Dpair_acc is None else Dpair_acc + Dp
        for v in vs:
            a = dA[v] / (np.linalg.norm(dA[v]) + 1e-12)
            b = dB[v] / (np.linalg.norm(dB[v]) + 1e-12)
            dself_acc.setdefault(v, []).append(1 - float(a @ b))
    Dpair = Dpair_acc / NSEED
    np.fill_diagonal(Dpair, np.inf)
    Dself = np.array([np.mean(dself_acc[v]) for v in variants])
    nn = Dpair.min(axis=1)  # nearest-sibling distance per variant
    identifiable = nn > Dself
    n = len(variants)
    # pair-level resolvable
    iu, ju = np.triu_indices(n, 1)
    pair_D = Dpair[iu, ju]
    pair_noise = 0.5 * (Dself[iu] + Dself[ju])
    resolvable_pairs = pair_D > pair_noise
    rows.append(dict(
        gene=g, n_var=n,
        frac_identifiable=round(float(identifiable.mean()), 3),
        frac_pairs_resolvable=round(float(resolvable_pairs.mean()), 3),
        median_nn_dist=round(float(np.median(nn)), 3),
        median_Dself=round(float(np.median(Dself)), 3),
        median_nn_over_Dself=round(float(np.median(nn / (Dself + 1e-12))), 3),
    ))

df = pd.DataFrame(rows)
df.to_csv(OUT_DIR / "pairwise_resolvability.csv", index=False)
print(df.to_string(index=False))
print(f"\nsaved -> {OUT_DIR / 'pairwise_resolvability.csv'}")
