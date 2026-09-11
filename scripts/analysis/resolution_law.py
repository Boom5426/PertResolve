#!/usr/bin/env python3
"""Calibrate what a benchmark's resolution buys, and invert it into a depth requirement.

Reads the table written by ``resolution_sweep.py`` and does three things.

**Chooses the axis by evidence.** Four candidate summaries of the measurement are scored
against every outcome: the mean squared separation over pairs (``rho2``, what Methods
asserts), the median and geometric mean of the nearest-competitor separation, and the
fraction of perturbations whose nearest competitor clears the noise. The one that actually
orders the outcomes is used; the others are reported so the choice is visible.

**Checks which outcomes are monotone before fitting anything to them.** The split-half
reproducibility reference should rise with resolution, while the probability of recovering the full
ordering of graded predictors need not: at the floor no predictor is distinguishable, and
at saturation the best two are both essentially perfect and again indistinguishable, so the
curve can turn over. Fitting a monotone calibration to a non-monotone outcome would produce
a calculator that is confidently wrong at the top, so monotonicity is tested and reported
rather than assumed.

**Inverts the calibration into a depth.** The nearest-competitor separation is a property
of the perturbations and does not change with sequencing depth; the noise does, as
``eta2 ~ 1/m``, which Phase 1 confirmed empirically. So the resolution attainable at depth
``m`` scales from a pilot at depth ``m0`` as

    rho2_nn(m) = rho2_nn(m0) * m / m0

and the depth needed to reach a target is read off the calibrated curve by inversion. The
calculator is deliberately built on the median rather than on the above-noise fraction: a
fraction is a threshold crossing and does not scale, so it cannot be extrapolated this way.

Usage:
  python resolution_law.py --sweep SWEEP_CSV --out OUT_DIR
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import reject_repo_results

_ap = argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--sweep", required=True, help="resolution_sweep.csv from resolution_sweep.py")
_ap.add_argument("--out", required=True, metavar="OUT_DIR",
                 help="directory receiving resolution_law.json and resolution_law.csv")
_args = _ap.parse_args()

OUT_DIR = reject_repo_results(_args.out)
OUT_DIR.mkdir(parents=True, exist_ok=True)
df = pd.read_csv(_args.sweep)

AXES = ["rho2", "rho2_nn_median", "rho2_nn_geomean", "frac_above_noise"]
OUTCOMES = ["ceiling_pds", "p_correct", "p_winner"]


def spearman(x, y) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return float("nan")
    return float(np.corrcoef(np.argsort(np.argsort(x[ok])),
                             np.argsort(np.argsort(y[ok])))[0, 1])


print("=== which measurement summary orders each outcome? ===")
print(f"{'axis':20s} {'n usable':>9s} " + "".join(f"{o:>14s}" for o in OUTCOMES))
axis_scores = {}
for axis in AXES:
    usable = int((np.isfinite(df[axis]) & np.isfinite(df.ceiling_pds)).sum())
    scores = [spearman(df[axis], df[o]) for o in OUTCOMES]
    axis_scores[axis] = scores[0]
    print(f"{axis:20s} {usable:9d} " + "".join(f"{s:14.3f}" for s in scores))

AXIS = max(axis_scores, key=lambda k: (-1 if np.isnan(axis_scores[k]) else axis_scores[k]))
print(f"\nchosen axis: {AXIS}")

print("\n=== is each outcome monotone in that axis? ===")
print("(Spearman on the lower and the upper half of the axis, separately: an outcome that")
print(" rises then falls shows a positive coefficient below and a negative one above.)")
sub = df[np.isfinite(df[AXIS])].copy()
cut = sub[AXIS].median()
monotone, turn = {}, {}
for o in OUTCOMES:
    lo = spearman(sub[sub[AXIS] <= cut][AXIS], sub[sub[AXIS] <= cut][o])
    hi = spearman(sub[sub[AXIS] > cut][AXIS], sub[sub[AXIS] > cut][o])
    monotone[o] = bool(np.isfinite(hi) and hi > 0)
    flag = "monotone" if monotone[o] else "TURNS OVER, not usable for inversion"
    print(f"  {o:14s} lower half {lo:+.3f}   upper half {hi:+.3f}   {flag}")
    if not monotone[o]:
        # locate the turnover, since where it happens decides whether it matters in the
        # range the real datasets occupy
        q = sub.groupby(pd.qcut(sub[AXIS], 12, duplicates="drop"),
                        observed=True)[[AXIS, o]].median()
        peak = q[o].idxmax()
        turn[o] = float(q.loc[peak, AXIS])
        ceil_at_peak = float(sub[np.isclose(sub[AXIS], turn[o], rtol=0.5)].ceiling_pds.median())
        print(f"                 peak near {AXIS} = {turn[o]:.4f}, "
              f"around an empirical reference of {ceil_at_peak:.3f}")

TARGETS = [o for o in OUTCOMES if monotone[o]]
if not TARGETS:
    raise SystemExit("no outcome is monotone in the chosen axis; nothing can be inverted")
print(f"\ninvertible outcomes: {', '.join(TARGETS)}")


def calibrate(frame: pd.DataFrame, axis: str, outcome: str, n_bins: int = 24):
    """Monotone calibration by binning, with a running maximum to enforce monotonicity.

    No functional form is assumed. The claim the sweep supports is that the outcome rises
    with the axis, not that it follows any particular curve, and a binned monotone
    interpolation is exactly that claim and no more.
    """
    d = frame[[axis, outcome]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(d) < n_bins * 2:
        n_bins = max(4, len(d) // 4)
    d = d.sort_values(axis)
    edges = np.quantile(d[axis], np.linspace(0, 1, n_bins + 1))
    edges = np.unique(edges)
    ax_vals = d[axis].to_numpy()
    out_vals = d[outcome].to_numpy()
    idx = np.clip(np.searchsorted(edges, ax_vals, side="right") - 1, 0, len(edges) - 2)
    xs, ys = [], []
    for b in range(len(edges) - 1):
        sel = idx == b
        if sel.sum() >= 3:
            xs.append(float(np.median(ax_vals[sel])))
            ys.append(float(np.median(out_vals[sel])))
    ys = list(np.maximum.accumulate(ys))
    return np.asarray(xs), np.asarray(ys)


curves = {}
rows = []
for outcome in TARGETS:
    for n_lo, n_hi, label in [(0, 30, "n<=30"), (30, 70, "30<n<=70"), (70, 10 ** 6, "n>70")]:
        frame = df[(df.n_var > n_lo) & (df.n_var <= n_hi)]
        if len(frame) < 20:
            continue
        xs, ys = calibrate(frame, AXIS, outcome)
        curves[f"{outcome}|{label}"] = dict(axis=AXIS, x=xs.tolist(), y=ys.tolist(),
                                            n_points=int(len(frame)))
        for x, y in zip(xs, ys):
            rows.append(dict(outcome=outcome, n_bucket=label, axis=AXIS, x=x, y=y))
        print(f"  calibrated {outcome:12s} for {label:9s} from {len(frame):4d} points, "
              f"{len(xs)} knots, range {ys.min():.3f} to {ys.max():.3f}")


def required_axis_value(outcome: str, n_var: int, target: float) -> float:
    """Smallest calibrated axis value reaching ``target``, or NaN if the curve never does."""
    label = "n<=30" if n_var <= 30 else ("30<n<=70" if n_var <= 70 else "n>70")
    key = f"{outcome}|{label}"
    if key not in curves:
        return float("nan")
    xs, ys = np.asarray(curves[key]["x"]), np.asarray(curves[key]["y"])
    hit = np.where(ys >= target)[0]
    return float(xs[hit[0]]) if len(hit) else float("nan")


def required_depth(pilot_rho2_nn_median: float, pilot_depth: int, n_var: int,
                   outcome: str, target: float, *, pilot_lo: float | None = None) -> float:
    """Cells per group needed to reach ``target``, from a pilot at ``pilot_depth``.

    The nearest-competitor separation belongs to the perturbations and does not change with
    depth; the noise does, as ``1 / m``. The attainable resolution therefore scales linearly
    in depth from the pilot, and the depth is read off the calibrated curve by inversion.

    Returns NaN, rather than a number, in three cases where a number would be fabricated:
    the calibration never reaches the target at this candidate-set size; the pilot estimate
    is not positive; or, when ``pilot_lo`` is supplied, its lower confidence bound does not
    clear zero. The last case is the one that matters in practice. A dataset whose
    separation is estimated at 0.0005 with an interval covering zero will otherwise be told
    to profile a few thousand cells per perturbation, which reads as a plan and is really an
    extrapolation from noise.
    """
    need = required_axis_value(outcome, n_var, target)
    if not np.isfinite(need) or pilot_rho2_nn_median <= 0:
        return float("nan")
    if pilot_lo is not None and not (pilot_lo > 0):
        return float("nan")
    return float(pilot_depth * need / pilot_rho2_nn_median)


print("\n=== the four real datasets, read through the calibration ===")
real = df[df.source == "real"]
prescriptions = []
for _, r in real.iterrows():
    for outcome in TARGETS:
        for target in (0.7, 0.9):
            lo = r.get("rho2_nn_median_lo", None)
            need = required_depth(r.rho2_nn_median, int(r.depth_m), int(r.n_var),
                                  outcome, target, pilot_lo=lo)
            prescriptions.append(dict(gene=r.gene, pilot_depth=int(r.depth_m),
                                      n_var=int(r.n_var),
                                      rho2_nn_median=r.rho2_nn_median,
                                      rho2_nn_median_lo=lo,
                                      outcome=outcome, target=target,
                                      required_cells_per_group=need,
                                      refused=bool(not np.isfinite(need))))
pres = pd.DataFrame(prescriptions)
print(pres.to_string(index=False))

pd.DataFrame(rows).to_csv(OUT_DIR / "resolution_law.csv", index=False)
pres.to_csv(OUT_DIR / "resolution_law_prescriptions.csv", index=False)
with open(OUT_DIR / "resolution_law.json", "w") as fh:
    json.dump(dict(axis=AXIS,
                   axis_spearman_vs_ceiling=axis_scores,
                   monotone=monotone,
                   turnover=turn,
                   curves=curves), fh, indent=2)
print(f"\nsaved -> {OUT_DIR / 'resolution_law.csv'}, "
      f"{OUT_DIR / 'resolution_law_prescriptions.csv'}, {OUT_DIR / 'resolution_law.json'}")
