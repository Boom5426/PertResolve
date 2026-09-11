#!/usr/bin/env python3
"""Measurement-resolution scaling for the four allele datasets, with uncertainty.

Supersedes ``floor_law_fit.py``, which is kept unchanged because it reproduces the
committed ``results/canonical/floor_law.csv``. Two defects in that estimator made its
x-axis unusable on exactly the datasets this study is about. The estimator itself now
lives in ``pertresolve.resolution.scaling``, which documents the model and is validated
against data with a known signal; this script supplies the cells and the reporting.

**1. The noise term was inflated by the wild-type mean.** ``floor_law_fit.py`` estimated
the per-profile noise from two half-profiles that subtract two *independent* wild-type
half-means, while the between-variant term was computed among profiles that all share one
wild-type mean, so that mean cancels there exactly. The wild-type sampling noise therefore
entered the noise term twice and the signal term not at all. Measured on the real splits
that term is as large as the variant sampling noise itself (TP53 at 25 cells per group:
36.9 against a true 37.6), so the noise estimate came out roughly doubled and ``delta2``
landed near -74 rather than near 0.

**2. The estimate was clipped at zero.** ``delta2 = max(pair2 - 2 * eta2, 0.0)`` turned
every negative estimate into an exact zero, which with defect 1 present was 13 of 16
gene-by-depth points. On data with no true signal, half of all unbiased estimates are
negative, so clipping reports certainty the data does not carry.

Corrected here, together with three additions the superseded script had no equivalent of:

* 95% bootstrap intervals over variants, the resampling unit used throughout this study;
* a permutation P value for the null that no between-variant signal exists. A
  likelihood-based variance-components fit (REML) was considered and rejected: it
  constrains the between variance to be non-negative, reintroducing exactly the boundary
  artefact that defect 2 is about, so it cannot serve as an independent check;
* disjoint cells for the two axes. ``floor_law_fit.py`` computed the between-variant term
  and the split-half reference from the same evaluation half, so one noise realisation entered
  both axes of the intended collapse plot. Each variant's cells are split here into four
  disjoint groups of ``m`` cells: two (Q, T) score the empirical reference and two (S1, S2) estimate
  the signal. Nothing is shared, at the cost of requiring 4m rather than 2m cells per
  variant, which is reported as ``n_var`` rather than hidden.

**Depth series on a fixed cohort.** Raising ``m`` also removes variants that no longer
have enough cells, so a depth trend over all eligible variants confounds depth with a
changing cohort. Rows are emitted twice: ``cohort=eligible`` uses every variant with at
least 4m cells, and ``cohort=fixed`` restricts the whole series to the variants eligible
at every depth, which is the comparison that isolates depth.

Usage:
  python resolution_scaling.py --out OUT_DIR [--base BASE] [--depths 25,50,100,150,250,400]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# The package is not installed by default; make this repository importable so that
# the shared path helpers and the estimator can be used when the script is run directly.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import (
    add_harness_to_path,
    reject_repo_results,
    require_inputs,
    resolve_base,
)
from pertresolve.resolution import (
    bootstrap_delta2,
    gram_of,
    permutation_null_delta2,
    signal_noise,
    stats_from_gram,
    tie_aware_pds,
)

_ap = argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--base", default=None,
                 help="directory holding the shared scorer and the gene "
                      "arrays (env: PERTRESOLVE_DATA)")
_ap.add_argument("--out", required=True, metavar="OUT_DIR",
                 help="directory receiving floor_law_v2.csv; may not be inside the "
                      "repository's results/")
_ap.add_argument("--depths", default="25,50,100,150,250,400",
                 help="comma-separated cells per group (default: %(default)s)")
_ap.add_argument("--n-seed", type=int, default=20,
                 help="independent cell splits averaged per point (default: %(default)s)")
_ap.add_argument("--n-boot", type=int, default=2000,
                 help="bootstrap resamples over variants (default: %(default)s)")
_ap.add_argument("--n-perm", type=int, default=2000,
                 help="permutations for the no-signal null (default: %(default)s)")
_ap.add_argument("--min-var", type=int, default=5,
                 help="minimum variants required to emit a row (default: %(default)s)")
_args = _ap.parse_args()

BASE = resolve_base(_args.base)
OUT_DIR = reject_repo_results(_args.out)
OUT_DIR.mkdir(parents=True, exist_ok=True)
require_inputs(BASE / "pertresolve_bench.csv")
add_harness_to_path(BASE)

import harness as H  # noqa: E402

DEPTHS = [int(x) for x in _args.depths.split(",")]
NSEED, NBOOT, NPERM, MIN_VAR = _args.n_seed, _args.n_boot, _args.n_perm, _args.min_var
H_SIGNAL = 2          # measurements used to estimate the signal; Q and T are extra
N_GROUPS = 2 + H_SIGNAL
PERM_SPLITS = 5       # splits the permutation null averages over, of the NSEED drawn
WT_TAGS = ("WT", "wt", "WT_control")

_BENCH = pd.read_csv(BASE / "pertresolve_bench.csv")
BENCH = {g: set(_BENCH[_BENCH.gene == g]["variant"]) for g in H.GENES}
GENE_CELLS = {g: (lambda X, lab: (X, np.asarray(lab)))(*H.load_gene(g)) for g in H.GENES}


def eligible(gene: str, m: int) -> list[str]:
    """Variants of ``gene`` in the curated benchmark with enough cells for N_GROUPS."""
    _, lab = GENE_CELLS[gene]
    return sorted(v for v in np.unique(lab)
                  if v not in WT_TAGS and v in BENCH[gene]
                  and int((lab == v).sum()) >= N_GROUPS * m)


def split_groups(gene: str, m: int, seed: int, variants: list[str]):
    """Q, T and the signal profiles of each variant, on mutually disjoint cells.

    Q and T each subtract their own wild-type half-mean, matching the convention the
    empirical reference is defined with. The signal profiles all subtract one shared wild-type
    mean, so it cancels in both the within-variant and the between-variant terms and
    cannot bias their difference.

    Returns Q, T of shape (n, K), S of shape (n, H_SIGNAL, K), the wild-type mean noise
    ``0.5 * ||w_Q - w_T||^2`` (the term the superseded estimator folded into its noise
    estimate, measured from the wild-type means themselves rather than from profile
    averages, which would be badly biased at JAK1's variant counts), and how many
    wild-type cells each reference mean was built from.

    The reference is always built from whatever wild-type cells exist rather than
    silently dropping to no subtraction at all: TP53 carries fewer than 200 wild-type
    cells, so at 100 cells per group a rule demanding 2m of them would leave the profiles
    as raw means, and a cosine score does not ignore the additive shift that removes.
    ``wt_per_ref`` records the depth actually used so a row that rests on a shallower
    reference is visible in the table.
    """
    X, lab = GENE_CELLS[gene]
    rng = np.random.RandomState(1000 + 7 * seed + H._GENE_SEED[gene])
    wt = np.where(np.isin(lab, WT_TAGS))[0]
    if len(wt) >= 2:
        rng.shuffle(wt)
        # Q and T need disjoint references to match the split-half-reference convention; the
        # signal groups share one, so it cancels out of both terms of delta2.
        per_ref = min(m, len(wt) // 2)
        wm_q = X[wt[:per_ref]].mean(0)
        wm_t = X[wt[per_ref:2 * per_ref]].mean(0)
        wm_s = wm_t
    else:
        per_ref = 0
        wm_q = wm_t = wm_s = np.zeros(X.shape[1], np.float32)

    Q, T, S = [], [], []
    for v in variants:
        idx = np.where(lab == v)[0]
        rng.shuffle(idx)
        Q.append(X[idx[:m]].mean(0) - wm_q)
        T.append(X[idx[m:2 * m]].mean(0) - wm_t)
        S.append([X[idx[(2 + h) * m:(3 + h) * m]].mean(0) - wm_s for h in range(H_SIGNAL)])
    wt_noise = 0.5 * float(np.sum((np.asarray(wm_q, np.float64)
                                   - np.asarray(wm_t, np.float64)) ** 2))
    return (np.stack(Q).astype(np.float64),
            np.stack(T).astype(np.float64),
            np.stack([np.stack(s) for s in S]).astype(np.float64),
            wt_noise, per_ref)


def run(gene: str, m: int, variants: list[str], cohort: str) -> dict | None:
    """One (gene, depth, cohort) row: signal, noise, their ratio, and the reference."""
    n = len(variants)
    if n < MIN_VAR:
        return None
    ssw_stack = np.zeros((NSEED, n))
    d2_stack = np.zeros((NSEED, n, n))
    ceil_acc, wt_acc, ref_acc, grams = [], [], [], []
    for seed in range(NSEED):
        Q, T, S, wt_noise, per_ref = split_groups(gene, m, seed, variants)
        ref_acc.append(per_ref)
        gram = gram_of(S)
        ssw, d2 = stats_from_gram(gram, n, H_SIGNAL)
        ssw_stack[seed] = ssw
        d2_stack[seed] = d2
        ceil_acc.append(tie_aware_pds(Q, T))
        wt_acc.append(wt_noise)
        if len(grams) < PERM_SPLITS:
            grams.append(gram)
    ssw_acc = ssw_stack.mean(axis=0)
    d2_acc = d2_stack.mean(axis=0)

    obs = signal_noise(ssw_acc, d2_acc, H_SIGNAL)
    eta2, delta2, rho2 = obs["eta2"], obs["delta2"], obs["rho2"]

    # Two levels: which variants were assayed, and which cells fell in which group.
    # Resampling variants alone gave intervals narrow enough to exclude zero for a
    # quantity that cannot be negative, which is a property of the interval, not the data.
    per_split = np.array([signal_noise(ssw_stack[s], d2_stack[s], H_SIGNAL)["delta2"]
                          for s in range(NSEED)])
    boot = bootstrap_delta2(ssw_stack, d2_stack, H_SIGNAL, n_boot=NBOOT, seed=0)
    d_lo, d_hi = np.percentile(boot, [2.5, 97.5])
    r_lo, r_hi = (d_lo / (2 * eta2), d_hi / (2 * eta2)) if eta2 > 0 else (np.nan, np.nan)

    null = permutation_null_delta2(grams, n, H_SIGNAL, n_perm=NPERM, seed=1)
    p_value = float((np.sum(null >= delta2) + 1) / (len(null) + 1))

    return dict(
        gene=gene, depth_m=m, cohort=cohort, n_var=n,
        eta2=round(eta2, 3),
        eta2_wt_excluded=round(float(np.mean(wt_acc)), 3),
        wt_cells_per_ref=int(np.min(ref_acc)),
        wt_ref_full_depth=bool(int(np.min(ref_acc)) >= m),
        perm_splits=len(grams),
        delta2=round(delta2, 3), delta2_lo=round(float(d_lo), 3), delta2_hi=round(float(d_hi), 3),
        delta2_sd_over_splits=round(float(per_split.std(ddof=1)), 3),
        rho2=round(rho2, 4), rho2_lo=round(float(r_lo), 4), rho2_hi=round(float(r_hi), 4),
        perm_p=round(p_value, 4),
        null_delta2_mean=round(float(null.mean()), 3),
        ceiling_pds=round(float(np.mean(ceil_acc)), 3),
    )


rows = []
for gene in H.GENES:
    per_depth = {m: eligible(gene, m) for m in DEPTHS}
    for m in DEPTHS:
        row = run(gene, m, per_depth[m], "eligible")
        if row:
            rows.append(row)
    usable = [set(v) for m, v in per_depth.items() if len(v) >= MIN_VAR]
    common = sorted(set.intersection(*usable)) if usable else []
    if len(common) >= MIN_VAR:
        for m in DEPTHS:
            if len(per_depth[m]) >= MIN_VAR:
                row = run(gene, m, common, "fixed")
                if row:
                    rows.append(row)
    else:
        print(f"{gene}: only {len(common)} variants are eligible at every usable depth, "
              f"below the minimum of {MIN_VAR}; eligible-cohort rows only")

df = pd.DataFrame(rows)
df.to_csv(OUT_DIR / "floor_law_v2.csv", index=False)
print(df.to_string(index=False))
print(f"\nsaved -> {OUT_DIR / 'floor_law_v2.csv'}")

elig = df[df.cohort == "eligible"]
sig = elig[elig.perm_p < 0.05]
print(f"\nBetween-variant signal above the permutation null (P < 0.05): "
      f"{len(sig)} of {len(elig)} eligible-cohort points")
if len(sig):
    print(sig[["gene", "depth_m", "n_var", "rho2", "rho2_lo", "rho2_hi", "perm_p",
               "ceiling_pds"]].to_string(index=False))
