"""Regenerate `results/pearson_delta_bootstrap_summary.csv` for Fig. 2c.

Fig. 2c plots a per-method direction-recovery interval. Rendering it with a
stochastic resample would make the panel non-deterministic, so the panel reads a
committed summary table instead. That is only defensible if the table itself is
reproducible, which is what this script provides: the summary shipped with the
figure refinement had no generator, and a plotted interval whose provenance is a
loose CSV is exactly the kind of claim that cannot be defended in review.

Recipe, verified to reproduce the shipped table exactly (all 19 methods, both
bounds, to floating-point equality):

  * source            results/results_v4_exttheta.csv
  * unit              one score row (method x gene x split x variant)
  * filter            finite `pearson_delta`; n is 518 for theta-only heads and
                      511 for the heads that also need an ESM embedding
  * statistic         mean of `pearson_delta` over the resampled rows
  * resamples         2,000
  * interval          percentile, 2.5 and 97.5
  * generator         one numpy default_rng(0) shared across methods, drawn in
                      the method order below, so the table is a single
                      deterministic stream rather than 19 independent ones

Caveat, stated because it affects interpretation and is not visible in the
table: the resampling unit is the score row, not the variant cluster. Rows for
the same variant under different splits are not independent, so these intervals
are anti-conservative if that dependence is strong. A cluster-aware sensitivity
analysis is the right follow-up before submission; this script deliberately
reproduces the shipped definition rather than silently substituting a different
one.

Usage:
    python scripts/figures/make_pearson_bootstrap_summary.py [--check]

``--check`` recomputes and compares against the committed table without writing,
exiting non-zero on any disagreement.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
SOURCE = REPO / "results" / "results_v4_exttheta.csv"
TARGET = REPO / "results" / "pearson_delta_bootstrap_summary.csv"

SEED = 0
N_RESAMPLES = 2000
PERCENTILES = (2.5, 97.5)

# Fixed order. The generator is shared across methods, so the order is part of
# the recipe: permuting it changes every interval.
METHODS = [
    "Ridge-theta", "Ridge-esm", "Ridge-esm+theta",
    "Lasso-theta", "Lasso-esm", "Lasso-esm+theta",
    "RF-theta", "RF-esm", "RF-esm+theta",
    "GBoost-theta", "GBoost-esm", "GBoost-esm+theta",
    "KNN-theta", "KNN-esm", "KNN-esm+theta",
    "MLP-theta", "MLP-esm", "MLP-esm+theta",
    "Gene-mean",
]


def compute() -> list[dict[str, object]]:
    """Return one summary row per method, in METHODS order."""
    frame = pd.read_csv(SOURCE)
    missing = [m for m in METHODS if not (frame.method == m).any()]
    if missing:
        raise SystemExit(f"{SOURCE.name} has no rows for: {', '.join(missing)}")

    rng = np.random.default_rng(SEED)
    rows = []
    for method in METHODS:
        values = frame.loc[frame.method == method, "pearson_delta"].dropna().to_numpy(float)
        if values.size == 0:
            raise SystemExit(f"{method}: no finite pearson_delta rows")
        index = rng.integers(0, values.size, size=(N_RESAMPLES, values.size))
        means = values[index].mean(axis=1)
        lo, hi = np.percentile(means, PERCENTILES)
        rows.append({
            "method": method,
            "pearson_delta": repr(float(values.mean())),
            "ci_lo": repr(float(lo)),
            "ci_hi": repr(float(hi)),
            "n": int(values.size),
        })
    return rows


def write(rows: list[dict[str, object]], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["method", "pearson_delta",
                                                    "ci_lo", "ci_hi", "n"])
        writer.writeheader()
        writer.writerows(rows)


def check(rows: list[dict[str, object]], path: Path) -> int:
    if not path.exists():
        print(f"FAIL {path} does not exist; run without --check to create it")
        return 1
    committed = {r["method"]: r for r in csv.DictReader(path.open(encoding="utf-8"))}
    bad = 0
    for row in rows:
        have = committed.get(row["method"])
        if have is None:
            print(f"FAIL {row['method']}: absent from the committed table")
            bad += 1
            continue
        for key in ("pearson_delta", "ci_lo", "ci_hi"):
            if abs(float(have[key]) - float(row[key])) > 1e-12:
                print(f"FAIL {row['method']}.{key}: committed {have[key]} "
                      f"recomputed {row[key]}")
                bad += 1
        if int(have["n"]) != int(row["n"]):
            print(f"FAIL {row['method']}.n: committed {have['n']} recomputed {row['n']}")
            bad += 1
    if bad == 0:
        print(f"PASS {path.name}: {len(rows)} methods reproduce exactly "
              f"(seed {SEED}, {N_RESAMPLES} resamples, percentile "
              f"{PERCENTILES[0]}/{PERCENTILES[1]})")
    return 1 if bad else 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="verify the committed table instead of writing it")
    args = parser.parse_args()
    rows = compute()
    if args.check:
        raise SystemExit(check(rows, TARGET))
    write(rows, TARGET)
    print(f"wrote {TARGET.relative_to(REPO)} ({len(rows)} methods)")


if __name__ == "__main__":
    main()
