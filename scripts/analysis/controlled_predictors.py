#!/usr/bin/env python3
"""1C: Controlled-predictor benchmark-resolution curve (KEYSTONE).

Synthetic predictors with a KNOWN true quality order, indexed by alpha in [0, 1]:

    pred(v, alpha) = (1-alpha) * gene_mean_build + alpha * build_delta[v]

build_delta[v] = pseudobulk delta of v from a BUILD subsample; the evaluation truth
eval_delta[v] is an independent DISJOINT subsample. alpha=0 is the gene-mean (no allele
information); alpha=1 is fully allele-specific (its ceiling is the finite-sample oracle,
not 1.0). True quality is monotone in alpha, so a valid benchmark must rank
PDS(alpha) monotone in alpha.

Benchmark RESOLUTION is not whether the point-estimate PDS(alpha) is monotone (with many
variants it always is), but whether the ordering SURVIVES the finite-test-set
uncertainty. We therefore bootstrap the held-out variant set (B resamples) and report:
  P_correct_order = P(no inversions in PDS(alpha) vs alpha)      [full-order recovery]
  P_winner        = P(argmax_alpha PDS == alpha=1)              [trust the leaderboard top]
  mean_tau        = mean Kendall tau(alpha, PDS)
against the measurement window (ceiling = PDS(alpha=1), the oracle at that depth).
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
                 help="directory receiving controlled_recovery.csv; may not be inside the "
                      "repository's results/")
_args = _ap.parse_args()

BASE = resolve_base(_args.base)
OUT_DIR = reject_repo_results(_args.out)
OUT_DIR.mkdir(parents=True, exist_ok=True)
add_harness_to_path(BASE)

import harness as H

NSEED = 10
NBOOT = 1000
WT_TAGS = ('WT', 'wt', 'WT_control')
ALPHAS = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
NA = len(ALPHAS)
NPAIR = NA * (NA - 1) // 2
DEPTHS = [25, 50, 100, 150, 250]
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


def inversions(yy):
    c = 0
    for i in range(NA):
        for j in range(i + 1, NA):
            if yy[i] >= yy[j]:
                c += 1
    return c


rows = []
for g in H.GENES:
    for m in DEPTHS:
        # per-variant, seed-averaged PDS for each alpha: pv[a] -> array over variants
        vs = None
        per_seed = {a: [] for a in ALPHAS}  # list of per-variant arrays
        for seed in range(NSEED):
            bd, ed = build_eval(g, m, seed)
            if bd is None or len(bd) < 5:
                continue
            cur = sorted(bd.keys())
            if vs is None:
                vs = cur
            if cur != vs:
                continue  # variant set stable across seeds by construction; skip if not
            Build = np.stack([bd[v] for v in vs])
            Eval = np.stack([ed[v] for v in vs])
            gm = Build.mean(0)
            Tn = norm(Eval)
            n = len(vs)
            for a in ALPHAS:
                Pred = (1 - a) * gm + a * Build
                D = 1 - norm(Pred) @ Tn.T
                per_seed[a].append(np.array([pds_row(D[i], i, n) for i in range(n)]))
        if vs is None or len(per_seed[1.0]) == 0:
            continue
        pv = {a: np.mean(per_seed[a], axis=0) for a in ALPHAS}  # seed-averaged per variant
        n = len(vs)
        ceiling = float(pv[1.0].mean())
        pds_point = {a: float(pv[a].mean()) for a in ALPHAS}
        # bootstrap over variants
        rng = np.random.RandomState(0)
        n_correct = n_winner = 0
        tau_sum = 0.0
        for _ in range(NBOOT):
            bi = rng.randint(0, n, n)
            yy = [pv[a][bi].mean() for a in ALPHAS]
            inv = inversions(yy)
            tau_sum += 1 - 2 * inv / NPAIR
            if inv == 0:
                n_correct += 1
            if int(np.argmax(yy)) == NA - 1:
                n_winner += 1
        rows.append(dict(
            gene=g, depth_m=m, n_var=n,
            ceiling_pds=round(ceiling, 3),
            pds_a0=round(pds_point[0.0], 3),
            pds_a1=round(pds_point[1.0], 3),
            P_correct_order=round(n_correct / NBOOT, 3),
            P_winner=round(n_winner / NBOOT, 3),
            mean_tau=round(tau_sum / NBOOT, 3),
        ))

df = pd.DataFrame(rows)
df.to_csv(OUT_DIR / "controlled_recovery.csv", index=False)
print(df.to_string(index=False))
print(f"\nsaved -> {OUT_DIR / 'controlled_recovery.csv'}")

print("\n=== resolution curve: measurement window (ceiling) vs ranking recovery ===")
d2 = df.copy()
d2['ceil_bin'] = pd.cut(d2['ceiling_pds'], [0, 0.52, 0.56, 0.65, 0.80, 1.01],
                        labels=['<=0.52', '0.52-0.56', '0.56-0.65', '0.65-0.80', '>0.80'])
print(d2.groupby('ceil_bin', observed=True).agg(
    n=('P_correct_order', 'size'),
    mean_P_correct=('P_correct_order', 'mean'),
    mean_P_winner=('P_winner', 'mean')).round(3).to_string())
