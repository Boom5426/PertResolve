#!/usr/bin/env python3
"""Trace benchmark resolution against the measurement, on real and dialled configurations.

Methods asserts that discrimination depends on the data only through the dimensionless
group ``rho`` and the competitor count. That assertion is false: at ``rho2`` held within
1% of 0.59, the attainable ceiling ranges from 0.999 to 0.632 across configurations that
differ only in the rank of their separation directions and the spread of their amplitudes
(``docs/RESULT_COLLAPSE_REFUTED_2026-08-04.md``).

The reason is that ``rho2`` is a mean over all pairs while a discrimination score is a
ranking statistic: a perturbation is identified when its own measurement is nearer than
every competitor, which only its **closest** competitor can spoil. Far-apart pairs inflate
the mean without ever being asked to be resolved. This script traces both quantities
against the ceiling so the substitution can be checked rather than asserted.

Points come from two sources and are labelled as such.

**Real** rows are the four allele datasets at several depths, measured with no
construction at all.

**Dialled** rows add a synthetic separation of controlled size and controlled geometry to
TP53 and KRAS cells, whose own between-variant signal is measured at ``|rho2| < 0.005``
with every interval covering zero (``results/canonical/floor_law_v2.csv``). Rescaling
measured profiles would not work as a dial, since it scales signal and noise together and
leaves their ratio fixed; the dial has to act on the signal alone. The directions come from
a fifth cell block that no scoring or estimation touches, cut from the same shuffle as the
scoring groups so disjointness is guaranteed. A direction correlated with the noise it is
later added to would produce the relationship by construction rather than measuring it.

Every variant's cells are split into five disjoint groups of ``m`` cells:

    D            direction block, used only when dialling
    Q, T         query and truth for the ceiling and for the graded predictors
    S1, S2       estimate the noise, the mean separation and the nearest-competitor one

The two axes of every point therefore come from different cells: the ceiling from Q and T,
the signal statistics from S1 and S2.

Nearest-competitor separations are formed per split and averaged across splits, then
summarised. Forming them from split-averaged profiles is wrong and silently catastrophic:
averaged profiles carry a fraction of the noise that the debiasing term describes, so every
separation is driven negative.

Two summaries of those separations are reported and one is usually undefined. The median
tolerates separations that debias below zero, which is the normal case near the floor, and
so does the fraction lying above the noise. A geometric mean does not exist once any
separation is non-positive and is reported as undefined rather than floored, since a
floored value is set by the floor rather than by the data while still printing as a
number.

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
    nearest_neighbour_rho2,
    signal_noise,
    stats_from_gram,
    summarise_separations,
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
                 help="genes supplying cells for the dialled rows (default: %(default)s)")
_ap.add_argument("--depths", default="25,50,100",
                 help="comma-separated cells per group (default: %(default)s)")
_ap.add_argument("--n-var", default="20,45,90",
                 help="candidate-set sizes for the dialled rows (default: %(default)s)")
_ap.add_argument("--rho2-targets", default="0,0.05,0.2,0.6,2.0,6.0",
                 help="mean signal-to-noise levels to dial (default: %(default)s)")
_ap.add_argument("--ranks", default="none,5,3,2,1",
                 help="rank of the dialled direction set; 'none' keeps it full "
                      "(default: %(default)s)")
_ap.add_argument("--heterogeneity", default="0,1,2",
                 help="lognormal s.d. of the dialled amplitudes (default: %(default)s)")
_ap.add_argument("--n-seed", type=int, default=10,
                 help="cell splits averaged per point (default: %(default)s)")
_ap.add_argument("--n-boot", type=int, default=500,
                 help="bootstraps for the ordering-recovery probability (default: %(default)s)")
_ap.add_argument("--pilot", action="store_true",
                 help="a small grid, to time the full sweep before committing to it")
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
RANKS = [None if x.strip().lower() == "none" else int(x) for x in _args.ranks.split(",")]
HETERO = [float(x) for x in _args.heterogeneity.split(",")]
NSEED, NBOOT = _args.n_seed, _args.n_boot
if _args.pilot:
    DEPTHS, N_VARS = DEPTHS[:1], N_VARS[:1]
    RHO2_TARGETS, RANKS, HETERO = [0.0, 0.6], [None, 2], [0.0, 2.0]
    NSEED, NBOOT = 3, 100

H_SIGNAL = 2
ALPHAS = [0.0, 0.5, 0.7, 0.85, 0.925, 0.96, 1.0]
WT_TAGS = ("WT", "wt", "WT_control")

_BENCH = pd.read_csv(BASE / "allele_perturb_bench.csv")
BENCH = {g: set(_BENCH[_BENCH.gene == g]["variant"]) for g in H.GENES}
GENE_CELLS = {g: (lambda X, lab: (X, np.asarray(lab)))(*H.load_gene(g)) for g in H.GENES}


def eligible(gene: str, m: int, n_groups: int) -> list[str]:
    _, lab = GENE_CELLS[gene]
    return sorted(v for v in np.unique(lab)
                  if v not in WT_TAGS and v in BENCH[gene]
                  and int((lab == v).sum()) >= n_groups * m)


def split_groups(gene: str, m: int, seed: int, variants: list[str], *, reserve_direction: bool):
    """Scoring groups, and the direction block when one is reserved, from one shuffle.

    Q and T each subtract their own wild-type half-mean, matching the convention the
    replicate ceiling is defined with. The signal groups share one, so it cancels from both
    terms of the signal estimate and cannot bias their difference.
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

    off = 1 if reserve_direction else 0
    D, Q, T, S = [], [], [], []
    for v in variants:
        idx = np.where(lab == v)[0]
        rng.shuffle(idx)
        if reserve_direction:
            D.append(X[idx[:m]].mean(0))
        Q.append(X[idx[off * m:(off + 1) * m]].mean(0) - wm_q)
        T.append(X[idx[(off + 1) * m:(off + 2) * m]].mean(0) - wm_t)
        S.append([X[idx[(off + 2 + h) * m:(off + 3 + h) * m]].mean(0) - wm_s
                  for h in range(H_SIGNAL)])

    U = None
    if reserve_direction:
        U = np.stack(D).astype(np.float64)
        U -= U.mean(axis=0, keepdims=True)      # a shift shared by all variants moves no score
        U /= np.linalg.norm(U, axis=1, keepdims=True) + 1e-12
    return (np.stack(Q).astype(np.float64),
            np.stack(T).astype(np.float64),
            np.stack([np.stack(s) for s in S]).astype(np.float64),
            U)


def reshape_directions(U: np.ndarray, rank, hetero: float, rng: np.random.RandomState):
    """Impose a rank and an amplitude spread, then restore the mean pairwise separation.

    Renormalising the mean is what makes the comparison informative: configurations differ
    in geometry at a fixed mean separation, so any change in the ceiling cannot be a change
    in the mean.
    """
    if rank is not None and rank < U.shape[0]:
        u, s, vt = np.linalg.svd(U, full_matrices=False)
        U = (u[:, :rank] * s[:rank]) @ vt[:rank]
    if hetero:
        U = U * rng.lognormal(0.0, hetero, size=U.shape[0])[:, None]
    gram = U @ U.T
    d2 = np.diag(gram)[:, None] + np.diag(gram)[None, :] - 2 * gram
    spread = float(d2[np.triu_indices(U.shape[0], 1)].mean())
    return U / np.sqrt(max(spread, 1e-12)) * np.sqrt(2.0)


def effective_rank(config: np.ndarray) -> float:
    """Participation ratio of the configuration's eigenvalue spectrum."""
    centred = config - config.mean(axis=0, keepdims=True)
    ev = np.linalg.eigvalsh(centred @ centred.T)
    ev = ev[ev > 1e-12]
    return float(ev.sum() ** 2 / (ev ** 2).sum()) if len(ev) else float("nan")


def graded_scores(Q: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Per-alpha, per-perturbation discrimination score for predictors of known quality.

    Predictor ``alpha`` predicts ``(1 - alpha) * panel_mean + alpha * Q_v``, so its true
    quality rises with alpha by construction. Q and T come from disjoint cells, so the best
    predictor's ceiling is the replicate ceiling rather than a perfect score.
    """
    panel_mean = Q.mean(axis=0)
    tn = T / (np.linalg.norm(T, axis=1, keepdims=True) + 1e-12)
    out = []
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
        out.append(np.asarray(scores))
    return np.stack(out)


def recovery(per_alpha: np.ndarray, seed: int = 0) -> tuple[float, float]:
    """Bootstrap over perturbations: full-order recovery and correct-winner rates."""
    rng = np.random.RandomState(seed)
    n = per_alpha.shape[1]
    correct = winner = 0
    for _ in range(NBOOT):
        means = per_alpha[:, rng.randint(0, n, n)].mean(axis=1)
        if all(means[k + 1] > means[k] for k in range(len(ALPHAS) - 1)):
            correct += 1
        if int(np.argmax(means)) == len(ALPHAS) - 1:
            winner += 1
    return correct / NBOOT, winner / NBOOT


def measure(gene: str, m: int, variants: list[str], label: str, *,
            dial=None, extra: dict | None = None) -> dict:
    """One row: both signal statistics, the ceiling and the ordering-recovery rates."""
    n = len(variants)
    ssw_acc = np.zeros(n)
    d2_acc = np.zeros((n, n))
    ceil, alpha_acc, configs = [], [], []
    per_split = []
    for seed in range(NSEED):
        Q, T, S, U = split_groups(gene, m, seed, variants,
                                  reserve_direction=dial is not None)
        if dial is not None:
            V = dial(U, seed)
            Q, T, S = Q + V, T + V, S + V[:, None, :]
            configs.append(effective_rank(V))
        ssw, d2 = stats_from_gram(gram_of(S), n, H_SIGNAL)
        ssw_acc += ssw
        d2_acc += d2
        ceil.append(tie_aware_pds(Q, T))
        alpha_acc.append(graded_scores(Q, T))
        per_split.append(S)
    obs = signal_noise(ssw_acc / NSEED, d2_acc / NSEED, H_SIGNAL)

    # nearest-competitor separations are formed per split against that configuration's own
    # noise, then averaged; forming them from averaged profiles over-subtracts and floors
    seps = np.mean([nearest_neighbour_rho2(S, obs["eta2"])["separations"]
                    for S in per_split], axis=0)
    nn = summarise_separations(seps, obs["eta2"])

    # The calculator must not extrapolate a depth from a resolution estimate that is
    # consistent with zero: no finite depth rescues a separation that is not there, and a
    # confident number would be fabricated. Bootstrap the median over perturbations so the
    # gate is available downstream.
    rng_nn = np.random.RandomState(11)
    boot_med = np.array([np.median(seps[rng_nn.randint(0, n, n)]) for _ in range(2000)])
    nn_lo, nn_hi = np.percentile(boot_med, [2.5, 97.5]) / (2.0 * obs["eta2"])

    p_correct, p_winner = recovery(np.mean(alpha_acc, axis=0))
    row = dict(
        source=label, gene=gene, depth_m=m, n_var=n,
        eta2=round(obs["eta2"], 3),
        rho2=round(obs["rho2"], 5),
        rho2_nn_median=round(nn["rho2_nn_median"], 5),
        rho2_nn_median_lo=round(float(nn_lo), 5),
        rho2_nn_median_hi=round(float(nn_hi), 5),
        rho2_nn_geomean=round(nn["rho2_nn_geomean"], 6),
        frac_above_noise=round(nn["frac_above_noise"], 4),
        nn_nonpositive=nn["n_nonpositive"],
        ceiling_pds=round(float(np.mean(ceil)), 4),
        p_correct=round(p_correct, 3), p_winner=round(p_winner, 3),
        eff_rank=round(float(np.mean(configs)), 2) if configs else np.nan,
    )
    row.update(extra or {})
    return row


rows = []

print("=== real datasets, measured with no construction ===")
for gene in H.GENES:
    for m in DEPTHS:
        variants = eligible(gene, m, 4)
        if len(variants) < 5:
            print(f"  {gene} m={m}: {len(variants)} eligible variants, skipped")
            continue
        r = measure(gene, m, variants, "real",
                    extra=dict(rho2_target=np.nan, rank=np.nan, hetero=np.nan))
        rows.append(r)
        print(f"  {gene:6s} m={m:3d} n={r['n_var']:3d} rho2={r['rho2']:8.4f} "
              f"rho2_nn_med={r['rho2_nn_median']:9.5f} above={r['frac_above_noise']:.2f} "
              f"ceiling={r['ceiling_pds']:.3f} P_order={r['p_correct']:.3f}")

print("\n=== dialled configurations on a substrate with no signal of its own ===")
for gene in SUBSTRATES:
    for m in DEPTHS:
        pool = eligible(gene, m, 5)
        for n_var in N_VARS:
            if len(pool) < n_var:
                print(f"  {gene} m={m}: {len(pool)} eligible, below the requested {n_var}")
                continue
            variants = list(np.random.RandomState(7).permutation(pool)[:n_var])
            Q0, T0, S0, U0 = split_groups(gene, m, 0, variants, reserve_direction=True)
            eta2_0 = signal_noise(*stats_from_gram(gram_of(S0), n_var, H_SIGNAL),
                                  H_SIGNAL)["eta2"]
            for rank in RANKS:
                for hetero in HETERO:
                    for target in RHO2_TARGETS:
                        amp = float(np.sqrt(max(target, 0.0) * 2.0 * eta2_0 / 2.0))

                        def dial(U, seed, rank=rank, hetero=hetero, amp=amp):
                            shaped = reshape_directions(U, rank, hetero,
                                                        np.random.RandomState(4242))
                            return shaped * amp

                        r = measure(gene, m, variants, "dialled", dial=dial,
                                    extra=dict(rho2_target=target,
                                               rank=-1 if rank is None else rank,
                                               hetero=hetero))
                        rows.append(r)
                        print(f"  {gene:5s} m={m:3d} n={n_var:3d} rank={str(rank):>4s} "
                              f"h={hetero:.0f} t={target:<4} rho2={r['rho2']:7.4f} "
                              f"rho2_nn_med={r['rho2_nn_median']:9.5f} "
                              f"above={r['frac_above_noise']:.2f} "
                              f"ceiling={r['ceiling_pds']:.3f} P_order={r['p_correct']:.3f}")

df = pd.DataFrame(rows)
df.to_csv(OUT_DIR / "resolution_sweep.csv", index=False)
print(f"\nsaved -> {OUT_DIR / 'resolution_sweep.csv'}  ({len(df)} points)")


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return float("nan")
    rx = np.argsort(np.argsort(x[ok]))
    ry = np.argsort(np.argsort(y[ok]))
    return float(np.corrcoef(rx, ry)[0, 1])


cands = ["rho2", "rho2_nn_median", "rho2_nn_geomean", "frac_above_noise"]
print("\nSpearman with the ceiling:")
print(f"{'subset':12s} {'n':>4s} " + "".join(f"{c:>20s}" for c in cands))
for name, sub in (("real", df[df.source == "real"]),
                  ("dialled", df[df.source == "dialled"]),
                  ("combined", df)):
    print(f"{name:12s} {len(sub):4d} "
          + "".join(f"{spearman(sub[c], sub.ceiling_pds):20.3f}" for c in cands))
n_censored = int(df.rho2_nn_geomean.isna().sum())
print(f"\ngeometric mean undefined (some separation debiases below zero): "
      f"{n_censored} of {len(df)} points")
