#!/usr/bin/env python3
"""Trace benchmark resolution continuously against the measurement signal-to-noise.

The four allele datasets give four points on the axis that decides whether a benchmark can
tell perturbations apart. Four points cannot establish that the relationship is a single
function of that axis, which is what Methods asserts. This traces the relationship
continuously, on real cells, by adding a synthetic separation of controlled size to a
substrate whose own separation is measured to be zero.

Why not simply rescale the measured profiles
--------------------------------------------
Multiplying a measured profile by a constant scales the signal and the noise together, so
the ratio of the two, the only thing that matters here, does not move at all. The dial has
to act on the signal alone.

Construction
------------
TP53 and KRAS are the substrate. Their between-variant signal is measured at
``rho2`` between -0.002 and 0.004 with every interval covering zero
(``results/canonical/floor_law_v2.csv``), so they supply realistic cells, realistic
covariance and realistic depth scaling with essentially no separation of their own.

Each variant's cells are split into five disjoint groups of ``m`` cells:

    D            estimates that variant's synthetic direction
    Q, T         query and truth for the replicate ceiling and the graded predictors
    S1, S2       estimate the signal and the noise

A unit direction ``u_v`` is built from D alone, so it is independent of every group used
for scoring or estimation; a direction taken from the scoring cells would correlate with
their noise and the collapse would follow by construction rather than by measurement.
``s * u_v`` is then added to Q, T, S1 and S2, which separates the variants by a controlled
amount without touching the noise.

Nothing downstream assumes the resulting separation is what the construction intended.
``rho2`` is measured from S1 and S2 with the estimator validated in
``tests/test_resolution_scaling.py``, on cells disjoint from Q and T, so the two axes of
every plotted point come from different cells.

Axes traced
-----------
``rho2``          measured signal-to-noise, dialled through ``s`` and through depth ``m``
``n_var``         number of candidate perturbations, the second argument of the law
``ceiling_pds``   discrimination attainable by a second measurement of the same variant
``p_correct``     probability that the benchmark recovers the true ordering of graded
                  predictors, and ``p_winner`` that it picks the best one. The predictors
                  interpolate between the panel mean and each variant's own measured
                  effect, so their true ordering is fixed by construction.

Depth and amplitude are varied independently. If discrimination depends on the data only
through ``rho2``, points reached by raising depth and points reached by raising amplitude
must fall on the same curve at matched ``n_var``; if they separate, the assertion is wrong.
That comparison is the point of the sweep. It is made from the emitted table rather than
reported as a column, because it is a statement about the whole surface rather than about
any one row.

Usage:
  python resolution_sweep.py --out OUT_DIR [--base BASE] [--pilot]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from alleleperturb.paths import (
    add_harness_to_path,
    reject_repo_results,
    require_inputs,
    resolve_base,
)
from alleleperturb.resolution import (
    gram_of,
    signal_noise,
    stats_from_gram,
    tie_aware_pds,
)

_ap = argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--base", default=None,
                 help="VCCompass compute workspace (env: VCCOMPASS_BASE)")
_ap.add_argument("--out", required=True, metavar="OUT_DIR",
                 help="directory receiving resolution_sweep.csv; may not be inside "
                      "the repository's results/")
_ap.add_argument("--substrates", default="TP53,KRAS",
                 help="genes supplying cells; their own signal must be near zero "
                      "(default: %(default)s)")
_ap.add_argument("--depths", default="25,50,100",
                 help="comma-separated cells per group (default: %(default)s)")
_ap.add_argument("--n-var", default="20,45,90",
                 help="comma-separated candidate-set sizes (default: %(default)s)")
_ap.add_argument("--rho2-targets", default="0,0.01,0.03,0.08,0.2,0.5,1.0,2.0,4.0,8.0",
                 help="signal-to-noise levels to dial (default: %(default)s)")
_ap.add_argument("--n-seed", type=int, default=10,
                 help="cell splits averaged per point (default: %(default)s)")
_ap.add_argument("--n-boot", type=int, default=500,
                 help="bootstraps for the ordering-recovery probability (default: %(default)s)")
_ap.add_argument("--pilot", action="store_true",
                 help="run a small grid to time the full sweep before committing to it")
_args = _ap.parse_args()

BASE = resolve_base(_args.base)
OUT_DIR = reject_repo_results(_args.out)
OUT_DIR.mkdir(parents=True, exist_ok=True)
require_inputs(BASE / "allele_perturb_bench.csv")
add_harness_to_path(BASE)

import harness as H  # noqa: E402

SUBSTRATES = _args.substrates.split(",")
DEPTHS = [int(x) for x in _args.depths.split(",")]
N_VARS = [int(x) for x in _args.n_var.split(",")]
RHO2_TARGETS = [float(x) for x in _args.rho2_targets.split(",")]
NSEED, NBOOT = _args.n_seed, _args.n_boot
if _args.pilot:
    DEPTHS, N_VARS = DEPTHS[:1], N_VARS[:1]
    RHO2_TARGETS = [0.0, 0.5, 4.0]
    NSEED, NBOOT = 3, 100

H_SIGNAL = 2
N_GROUPS = 1 + 2 + H_SIGNAL          # D, Q, T, S1, S2
ALPHAS = [0.0, 0.5, 0.7, 0.85, 0.925, 0.96, 1.0]
WT_TAGS = ("WT", "wt", "WT_control")

_BENCH = pd.read_csv(BASE / "allele_perturb_bench.csv")
BENCH = {g: set(_BENCH[_BENCH.gene == g]["variant"]) for g in H.GENES}
GENE_CELLS = {g: (lambda X, lab: (X, np.asarray(lab)))(*H.load_gene(g)) for g in SUBSTRATES}


def eligible(gene: str, m: int) -> list[str]:
    _, lab = GENE_CELLS[gene]
    return sorted(v for v in np.unique(lab)
                  if v not in WT_TAGS and v in BENCH[gene]
                  and int((lab == v).sum()) >= N_GROUPS * m)


def split_groups(gene: str, m: int, seed: int, variants: list[str]):
    """Direction block and the four scoring groups, cut from one shuffle per variant.

    The five blocks come from a single permutation of each variant's cells, which is what
    makes them disjoint. Drawing the direction block from its own shuffle, as an earlier
    draft did, leaves it overlapping the scoring cells, and a direction correlated with
    the noise it is later added to would produce the collapse by construction instead of
    measuring it.

    The direction is re-estimated per split, from that split's own reserved block, and is
    centred on the panel mean: a shift shared by every variant moves no discrimination
    score, so only the part that separates variants is kept.

    Returns Q, T of shape (n, K), S of shape (n, H_SIGNAL, K) and unit directions U of
    shape (n, K).
    """
    X, lab = GENE_CELLS[gene]
    rng = np.random.RandomState(1000 + 7 * seed + H._GENE_SEED[gene])
    wt = np.where(np.isin(lab, WT_TAGS))[0]
    if len(wt) >= 2:
        rng.shuffle(wt)
        per_ref = min(m, len(wt) // 2)
        wm_q, wm_t = X[wt[:per_ref]].mean(0), X[wt[per_ref:2 * per_ref]].mean(0)
        wm_s = wm_t
    else:
        wm_q = wm_t = wm_s = np.zeros(X.shape[1], np.float32)

    D, Q, T, S = [], [], [], []
    for v in variants:
        idx = np.where(lab == v)[0]
        rng.shuffle(idx)
        D.append(X[idx[:m]].mean(0))
        Q.append(X[idx[m:2 * m]].mean(0) - wm_q)
        T.append(X[idx[2 * m:3 * m]].mean(0) - wm_t)
        S.append([X[idx[(3 + h) * m:(4 + h) * m]].mean(0) - wm_s for h in range(H_SIGNAL)])

    U = np.stack(D).astype(np.float64)
    U -= U.mean(axis=0, keepdims=True)
    U /= np.linalg.norm(U, axis=1, keepdims=True) + 1e-12
    return (np.stack(Q).astype(np.float64),
            np.stack(T).astype(np.float64),
            np.stack([np.stack(s) for s in S]).astype(np.float64),
            U)


def graded_recovery(Q: np.ndarray, T: np.ndarray):
    """How often the benchmark recovers the known ordering of graded predictors.

    Predictor ``alpha`` predicts ``(1 - alpha) * panel_mean + alpha * Q_v``, so its true
    quality rises with alpha by construction. ``Q`` and ``T`` come from disjoint cells, so
    the best predictor's ceiling is the replicate ceiling rather than a perfect score.
    """
    panel_mean = Q.mean(axis=0)
    tn = T / (np.linalg.norm(T, axis=1, keepdims=True) + 1e-12)
    per_alpha = []
    for a in ALPHAS:
        pred = (1 - a) * panel_mean + a * Q
        pn = pred / (np.linalg.norm(pred, axis=1, keepdims=True) + 1e-12)
        dist = 1.0 - pn @ tn.T
        n = dist.shape[1]
        scores = []
        for i in range(dist.shape[0]):
            row = dist[i]
            target = row[i]
            less = int((row < target - 1e-12).sum())
            eq = int((np.abs(row - target) <= 1e-12).sum())
            scores.append(1.0 - (less + (eq - 1) / 2.0) / (n - 1))
        per_alpha.append(np.asarray(scores))
    return np.stack(per_alpha)                     # (n_alpha, n_var)


def recovery_probabilities(per_alpha: np.ndarray, rng: np.random.RandomState):
    """Bootstrap over perturbations: full-order recovery and correct-winner rates."""
    n = per_alpha.shape[1]
    correct = winner = 0
    for _ in range(NBOOT):
        idx = rng.randint(0, n, n)
        means = per_alpha[:, idx].mean(axis=1)
        if all(means[k + 1] > means[k] for k in range(len(ALPHAS) - 1)):
            correct += 1
        if int(np.argmax(means)) == len(ALPHAS) - 1:
            winner += 1
    return correct / NBOOT, winner / NBOOT


rows = []
for gene in SUBSTRATES:
    for m in DEPTHS:
        pool = eligible(gene, m)
        for n_var in N_VARS:
            if len(pool) < n_var:
                print(f"{gene} m={m}: {len(pool)} variants eligible, "
                      f"below the requested {n_var}; skipped")
                continue
            variants = list(np.random.RandomState(7).permutation(pool)[:n_var])

            # noise level and direction spread, needed to convert a target rho2 into an
            # amplitude; measured once per configuration on the first split
            Q0, T0, S0, U0 = split_groups(gene, m, 0, variants)
            ssw0, d20 = stats_from_gram(gram_of(S0), n_var, H_SIGNAL)
            eta2_0 = signal_noise(ssw0, d20, H_SIGNAL)["eta2"]
            iu = np.triu_indices(n_var, 1)
            gu = U0 @ U0.T
            spread = float((np.diag(gu)[:, None] + np.diag(gu)[None, :] - 2 * gu)[iu].mean())

            for target in RHO2_TARGETS:
                amp = float(np.sqrt(max(target, 0.0) * 2.0 * eta2_0 / max(spread, 1e-12)))
                ssw_acc = np.zeros(n_var)
                d2_acc = np.zeros((n_var, n_var))
                ceil_acc, alpha_acc = [], []
                for seed in range(NSEED):
                    Q, T, S, U = split_groups(gene, m, seed, variants)
                    Uv = U * amp
                    Q, T = Q + Uv, T + Uv
                    S = S + Uv[:, None, :]
                    ssw, d2 = stats_from_gram(gram_of(S), n_var, H_SIGNAL)
                    ssw_acc += ssw
                    d2_acc += d2
                    ceil_acc.append(tie_aware_pds(Q, T))
                    alpha_acc.append(graded_recovery(Q, T))
                obs = signal_noise(ssw_acc / NSEED, d2_acc / NSEED, H_SIGNAL)
                p_correct, p_winner = recovery_probabilities(
                    np.mean(alpha_acc, axis=0), np.random.RandomState(0))
                rows.append(dict(
                    substrate=gene, depth_m=m, n_var=n_var,
                    rho2_target=target, amplitude=round(amp, 4),
                    eta2=round(obs["eta2"], 3),
                    rho2=round(obs["rho2"], 4),
                    ceiling_pds=round(float(np.mean(ceil_acc)), 4),
                    p_correct=round(p_correct, 3), p_winner=round(p_winner, 3),
                    pds_alpha0=round(float(np.mean(alpha_acc, axis=0)[0].mean()), 4),
                    pds_alpha1=round(float(np.mean(alpha_acc, axis=0)[-1].mean()), 4),
                ))
                print(f"  {gene} m={m:3d} n={n_var:3d} target={target:<5} "
                      f"rho2={obs['rho2']:.4f} ceiling={np.mean(ceil_acc):.3f} "
                      f"P_order={p_correct:.3f}")

df = pd.DataFrame(rows)
df.to_csv(OUT_DIR / "resolution_sweep.csv", index=False)
print(f"\nsaved -> {OUT_DIR / 'resolution_sweep.csv'}")
print(f"{len(df)} points, rho2 from {df.rho2.min():.4f} to {df.rho2.max():.3f}")
