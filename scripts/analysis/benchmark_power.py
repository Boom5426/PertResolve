#!/usr/bin/env python3
"""Arm D: how large must a true quality difference be before a benchmark of N variants
recovers it?

Protocol frozen in ``docs/PREREG_CANDIDATE_AND_POWER_v1_ARM_D.md`` before this ran. Nothing
here may be revised after a result is seen.

The question is deliberately not "how many variants does a benchmark need". That number does
not exist and the paper's own argument says so. It is: at a given measurement resolution and
a given benchmark size, how large a difference between two predictors is recoverable.

Two things make this different from the recovery probabilities already in Fig. 5d.

**The axis is the realised gap, not the interpolation weight.** The map from weight to score
saturates hard, by a median factor of 8.3 and up to 40 across the committed ladders
(``scripts/analysis/alpha_ladder_audit.py``), so a curve keyed on weights reports a different
quantity than it claims. Everything below is keyed on
``ΔPDS = mean(pv[a_hi]) - mean(pv[a_lo])`` at the full cohort.

**Subsets are drawn without replacement.** Fig. 5d bootstraps the existing cohort n-of-n,
which answers how stable that cohort's verdict is. Drawing a subset of size ``N_test``
answers what a benchmark of that size would have concluded, which is the design question.
The comparison is paired: both predictors are scored on the same subset.

Usage:
    python scripts/analysis/benchmark_power.py [--base BASE] --out OUT

    --base   directory holding the per-gene arrays (env: PERTRESOLVE_DATA).
    --out    directory receiving benchmark_power.csv; may not be inside results/.
"""
from __future__ import annotations

import argparse
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import (  # noqa: E402
    add_harness_to_path,
    reject_repo_results,
    require_inputs,
    resolve_base,
)

# Frozen constants. ALPHAS, DEPTHS and NSEED match controlled_predictors.py so this arm and
# Fig. 5d score the same design-ordered family.
NSEED = 10
ALPHAS = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
DEPTHS = (25, 50, 100, 150, 250)
N_TEST_LADDER = (5, 10, 20, 40, 80, 160)
N_DRAW = 200
TOL = 1e-12
WT_TAGS = ("WT", "wt", "WT_control")


def norm(M: np.ndarray) -> np.ndarray:
    return M / (np.linalg.norm(M, axis=1, keepdims=True) + 1e-12)


def pds_row(drow: np.ndarray, ci: int, n: int) -> float:
    """Canonical tie-aware mid-rank PDS."""
    td = drow[ci]
    if not np.isfinite(td):
        return 0.5
    less = int((drow < td - TOL).sum())
    eq = int((np.abs(drow - td) <= TOL).sum())
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1) if n > 1 else 0.5


def build_eval(H, gene_cells, g: str, m: int, seed: int):
    """Disjoint build and evaluation halves, verbatim from controlled_predictors.py.

    Each variant contributes ``m`` cells to the build profile and a disjoint ``m`` to the
    evaluation target, and the two are referenced to **separate** wild-type halves so the
    predictor and its target share no cells and no control draw.
    """
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


def per_variant_scores(H, gene_cells, g: str, m: int):
    """Seed-averaged per-variant PDS for every weight, or None if the condition is empty."""
    vs, per_seed = None, {a: [] for a in ALPHAS}
    for seed in range(NSEED):
        bd, ed = build_eval(H, gene_cells, g, m, seed)
        if bd is None or len(bd) < 5:
            continue
        cur = sorted(bd)
        if vs is None:
            vs = cur
        if cur != vs:
            continue
        Build = np.stack([bd[v] for v in vs])
        Tn = norm(np.stack([ed[v] for v in vs]))
        gm = Build.mean(0)
        n = len(vs)
        for a in ALPHAS:
            D = 1 - norm((1 - a) * gm + a * Build) @ Tn.T
            per_seed[a].append(np.array([pds_row(D[i], i, n) for i in range(n)]))
    if vs is None or not per_seed[1.0]:
        return None, None
    return vs, {a: np.mean(per_seed[a], axis=0) for a in ALPHAS}


def power_rows(gene: str, depth: int, vs, pv, rng) -> list[dict]:
    """One row per (weight pair, N_test): the realised gap and the recovery probability."""
    n = len(vs)
    ladder = [k for k in N_TEST_LADDER if k <= n] + [n]
    rows = []
    for a_lo, a_hi in combinations(ALPHAS, 2):
        lo, hi = pv[a_lo], pv[a_hi]
        delta = float(hi.mean() - lo.mean())          # realised ΔPDS on the full cohort
        for k in ladder:
            if k == n:
                # The whole cohort is one draw, and the gate below checks it.
                p = float(hi.mean() > lo.mean())
                n_draw = 1
            else:
                idx = np.array([rng.choice(n, k, replace=False) for _ in range(N_DRAW)])
                p = float((hi[idx].mean(axis=1) > lo[idx].mean(axis=1)).mean())
                n_draw = N_DRAW
            rows.append(dict(gene=gene, depth_m=depth, n_cohort=n, alpha_lo=a_lo,
                             alpha_hi=a_hi, delta_pds=round(delta, 6), n_test=k,
                             p_correct=round(p, 5), n_draw=n_draw,
                             pds_lo=round(float(lo.mean()), 5),
                             pds_hi=round(float(hi.mean()), 5)))
    return rows


def check_gate(df: pd.DataFrame) -> None:
    """At the full cohort, recovery must be the indicator that the realised gap is positive."""
    full = df[df.n_test == df.n_cohort]
    expect = (full.delta_pds > 0).astype(float)
    bad = int((full.p_correct.to_numpy() != expect.to_numpy()).sum())
    if bad:
        raise SystemExit(f"gate failed: {bad} full-cohort rows disagree with the sign of "
                         "their own realised gap; the pairing is wrong")
    print(f"PASS  full-cohort recovery equals sign(ΔPDS) on all {len(full)} rows")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=None)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    base = resolve_base(args.base)
    out_dir = reject_repo_results(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    require_inputs(base / "pertresolve_bench.csv")
    add_harness_to_path(base)

    import harness as H  # noqa: E402

    gene_cells = {g: (lambda X, lab: (X, np.asarray(lab)))(*H.load_gene(g))
                  for g in H.GENES}
    rng = np.random.default_rng(20260828)
    rows = []
    for g in H.GENES:
        for m in DEPTHS:
            vs, pv = per_variant_scores(H, gene_cells, g, m)
            if vs is None:
                print(f"  {g} m={m}: no evaluable cohort, skipped")
                continue
            rows += power_rows(g, m, vs, pv, rng)
            print(f"  {g} m={m}: n={len(vs)}, ceiling {pv[1.0].mean():.3f}", flush=True)

    df = pd.DataFrame(rows)
    check_gate(df)
    df.to_csv(out_dir / "benchmark_power.csv", index=False)
    print(f"\nwrote {out_dir / 'benchmark_power.csv'} ({len(df):,} rows)")


if __name__ == "__main__":
    main()
