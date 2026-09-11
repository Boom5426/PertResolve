#!/usr/bin/env python3
"""Top-k allele shortlisting from the canonical ranking, for printed Figure 2h.

Every convention here is frozen by docs/FIG2_TOPK_ANALYSIS_LOCK_v1.md and none of them
may be revisited after seeing an output. Read that file before changing anything in this
one; in particular section 0 records why this is a lock and not a pre-registration.

The question is a tail statistic of a ranking that already exists. PDS summarises the
same-gene ranking of a held-out allele by its tie-aware mid-rank; Hit@k asks instead
whether the true allele lands in the top k of that same list. Nothing is re-predicted,
re-featurised or re-scored: this script consumes the raw ranking counts that
``train_only_grid.py --dump-rank-counts`` writes, which are the three numbers
``pds_row`` itself reduces to a score.

    Hit@k = clip(k - less, 0, eq) / eq          (expected hit under random tie-breaking)
    null  = min(k / n, 1)                        (per query, using that query's own pool)

The tie-aware form is not cosmetic. A fully tied predictor has less = 0 and eq = n, so it
scores exactly k/n, which is its own null. Under a hard ``rank <= k`` readout the same
predictor scores zero and would be drawn below chance, which would be a statement about
the readout rather than about the predictor.

Aggregation follows the canonical PDS estimator exactly: mean over the 15 cell-subsample
seeds, then over splits, then over distinct (gene, variant) queries, with a percentile
bootstrap over those queries. One bootstrap draw scores all four k, so the recovery curve
keeps its covariance.

Usage:
    python3 scripts/analysis/fig2_topk.py \\
        --counts <dir>/per_seed_rank_counts_trainonly.csv.gz \\
        --pds    results/canonical/per_seed_variant_pds.csv.gz \\
        --out    <scratch dir>
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import reject_repo_results  # noqa: E402

#: Frozen by the lock, section 2.3. The measurement-arm panel
#: results/canonical/split_half_topk.csv (printed Fig. 3f) uses the same four values, and
#: the two arms are meant to be read against each other.
K_VALUES = (1, 3, 5, 10)

#: Frozen by the lock, section 2.7. Matches train_only_grid.py.
NBOOT = 2000
BOOT_SEED = 0

#: Frozen by the lock, section 2.2. Asserted, not assumed.
POOL_SIZES = {"TP53": 98, "KRAS": 92, "GATA1": 254, "JAK1": 26}

TOL = 1e-12


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True,
                    help="per_seed_rank_counts_trainonly.csv.gz from train_only_grid.py "
                         "--dump-rank-counts")
    ap.add_argument("--pds", required=True,
                    help="the committed per_seed_variant_pds.csv.gz, used as the "
                         "acceptance gate of the lock's section 2.1")
    ap.add_argument("--out", required=True, metavar="OUT_DIR",
                    help="directory receiving the tables; may not be inside results/")
    return ap.parse_args()


def hit_at_k(less: np.ndarray, eq: np.ndarray, k: int) -> np.ndarray:
    """Expected probability of landing in the top k when ties are broken at random."""
    return np.clip(k - less, 0, eq) / eq


def null_at_k(n: np.ndarray, k: int) -> np.ndarray:
    return np.minimum(k / n, 1.0)


def acceptance_gate(counts: pd.DataFrame, pds_path: Path) -> None:
    """The re-run must reproduce the committed scores exactly, or the analysis stops.

    This is the lock's section 2.1. It exists because the committed table is produced by
    refitting the heads on train-only rescaled arrays, and scoring the stored preds5
    predictions instead disagrees on 21.7% of rows. A mean tolerates that; a tail
    statistic does not.
    """
    committed = pd.read_csv(pds_path)
    key = ["method", "gene", "split", "variant", "seed"]
    merged = counts.merge(committed, on=key, how="inner", suffixes=("", "_committed"))
    if len(merged) != len(committed):
        raise SystemExit(
            f"acceptance gate: the re-run covers {len(merged)} of the committed "
            f"{len(committed)} rows; the ranking source is not the canonical one")
    diff = np.abs(merged["pds_recomputed"] - merged["pds"]).max()
    if diff > 5e-7:
        raise SystemExit(
            f"acceptance gate: max |pds difference| = {diff:.3e} against the committed "
            "table. The lock forbids adjusting the analysis to fit; find the divergence.")
    print(f"acceptance gate PASS: {len(merged)} rows, max |pds diff| = {diff:.3e} "
          "(committed table is printed to six decimals)")


def contract_tests(rows: pd.DataFrame) -> None:
    """The lock's section 3, asserted rather than inspected."""
    for gene, expected in POOL_SIZES.items():
        seen = sorted(rows.loc[rows.gene == gene, "n"].unique())
        if seen != [expected]:
            raise SystemExit(f"pool size for {gene} is {seen}, not [{expected}]")

    finite = rows[rows.finite]
    # 7. the dumped counts and the scored value come from one computation
    reconstructed = 1.0 - (finite["less"] + (finite["eq"] - 1) / 2.0) / (finite["n"] - 1)
    gap = np.abs(reconstructed - finite["pds_recomputed"]).max()
    if gap > 1e-12:
        raise SystemExit(f"counts do not reproduce their own pds: max gap {gap:.3e}")

    for k in K_VALUES:
        h = rows[f"hit_{k}"].to_numpy()
        if h.min() < -TOL or h.max() > 1 + TOL:
            raise SystemExit(f"hit_{k} leaves [0, 1]: {h.min()} to {h.max()}")
        expected_null = null_at_k(rows["n"].to_numpy(), k)
        if np.abs(rows[f"null_{k}"].to_numpy() - expected_null).max() > TOL:
            raise SystemExit(f"null_{k} is not min(k/n, 1) row for row")
    for a, b in zip(K_VALUES, K_VALUES[1:]):
        if (rows[f"hit_{a}"] > rows[f"hit_{b}"] + TOL).any():
            raise SystemExit(f"hit is not monotone in k between {a} and {b}")
        if (rows[f"null_{a}"] > rows[f"null_{b}"] + TOL).any():
            raise SystemExit(f"null is not monotone in k between {a} and {b}")

    # 6. a fully tied row sits exactly on its null
    tied = finite[finite["eq"] == finite["n"]]
    if len(tied):
        for k in K_VALUES:
            off = np.abs(tied[f"hit_{k}"] - null_at_k(tied["n"].to_numpy(), k)).max()
            if off > TOL:
                raise SystemExit(
                    f"a fully tied row is {off:.3e} off its null at k={k}; the tie-aware "
                    "convention is the whole reason a tied predictor must land on chance")
        print(f"contract: {len(tied)} fully tied rows sit exactly on their null")
    print("contract tests PASS")


def aggregate(rows: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Seeds, then splits, then distinct (gene, variant). The canonical order."""
    per_seed = rows.groupby(["method", "gene", "split", "variant"], observed=True)[columns].mean()
    per_query = per_seed.reset_index().groupby(
        ["method", "gene", "variant"], observed=True)[columns].mean().reset_index()
    return per_query


def bootstrap(per_query: pd.DataFrame, columns: list[str]) -> dict[str, tuple[float, float]]:
    """One draw scores every column, so the recovery curve keeps its covariance."""
    values = per_query[columns].to_numpy(float)
    rng = np.random.default_rng(BOOT_SEED)
    n = len(values)
    draws = np.empty((NBOOT, values.shape[1]))
    for b in range(NBOOT):
        draws[b] = values[rng.integers(0, n, n)].mean(axis=0)
    lo, hi = np.percentile(draws, [2.5, 97.5], axis=0)
    return {c: (float(lo[i]), float(hi[i])) for i, c in enumerate(columns)}


def main() -> None:
    args = parse_args()
    out_dir = reject_repo_results(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    counts = pd.read_csv(args.counts)
    counts = counts.rename(columns={"pds": "pds_recomputed"}) if "pds" in counts else counts

    # The counts file carries no score, so recompute it here from the counts alone. That
    # is deliberate: the gate below then compares two independent routes to the same
    # number rather than a value against itself.
    finite_mask = counts["finite"].astype(bool)
    counts["pds_recomputed"] = 0.5
    f = finite_mask.to_numpy()
    counts.loc[f, "pds_recomputed"] = 1.0 - (
        counts.loc[f, "less"] + (counts.loc[f, "eq"] - 1) / 2.0) / (counts.loc[f, "n"] - 1)
    counts["finite"] = f

    acceptance_gate(counts, Path(args.pds))

    less = counts["less"].to_numpy(float)
    eq = counts["eq"].to_numpy(float)
    n = counts["n"].to_numpy(float)
    for k in K_VALUES:
        null_k = null_at_k(n, k)
        counts[f"null_{k}"] = null_k
        hit = np.where(f, hit_at_k(less, np.maximum(eq, 1), k), null_k)
        # Lock 2.6: a prediction with no ranking is uninformative, so it takes its own
        # query's null expectation. Scoring it as a miss would turn "no ordering
        # information" into "negative ordering information".
        counts[f"hit_{k}"] = hit
        # Lock 2.6: the mandatory sensitivity arm, computed here so it can never be
        # produced later by a different route.
        counts[f"hit0_{k}"] = np.where(f, hit, 0.0)

    contract_tests(counts)

    hit_cols = [f"hit_{k}" for k in K_VALUES]
    hit0_cols = [f"hit0_{k}" for k in K_VALUES]
    null_cols = [f"null_{k}" for k in K_VALUES]
    per_query = aggregate(counts, hit_cols + hit0_cols + null_cols + ["n"])
    per_query["n"] = per_query["n"].round().astype(int)

    for k in K_VALUES:
        per_query[f"paired_{k}"] = per_query[f"hit_{k}"] - per_query[f"null_{k}"]
    paired_cols = [f"paired_{k}" for k in K_VALUES]

    coverage = counts.groupby("method", observed=True)["finite"].agg(
        n_rows="size", n_finite="sum")
    coverage["n_nonfinite"] = coverage["n_rows"] - coverage["n_finite"]

    summary = []
    for method, sub in per_query.groupby("method", observed=True):
        ci = bootstrap(sub, hit_cols + hit0_cols + paired_cols)
        cov = coverage.loc[method]
        for k in K_VALUES:
            lo, hi = ci[f"hit_{k}"]
            plo, phi = ci[f"paired_{k}"]
            summary.append(dict(
                method=method, k=k, n_queries=len(sub),
                hit_rate=float(sub[f"hit_{k}"].mean()), ci_lo=lo, ci_hi=hi,
                null_rate=float(sub[f"null_{k}"].mean()),
                paired_diff=float(sub[f"paired_{k}"].mean()),
                paired_lo=plo, paired_hi=phi,
                above_null=bool(plo > 0), below_null=bool(phi < 0),
                hit_rate_zero_as_failure=float(sub[f"hit0_{k}"].mean()),
                zero_ci_lo=ci[f"hit0_{k}"][0], zero_ci_hi=ci[f"hit0_{k}"][1],
                n_scored_rows=int(cov["n_rows"]), n_nonfinite_rows=int(cov["n_nonfinite"]),
            ))
    summary = pd.DataFrame(summary).sort_values(["method", "k"])

    ranking = counts[["method", "gene", "split", "variant", "seed",
                      "less", "eq", "n", "finite", "pds_recomputed"]].rename(
        columns={"pds_recomputed": "pds"})
    ranking.to_csv(out_dir / "fig2_ranking_per_query.csv.gz", index=False,
                   float_format="%.6f", compression="gzip")
    per_query.to_csv(out_dir / "fig2_topk_per_query.csv", index=False, float_format="%.6f")
    summary.to_csv(out_dir / "fig2_topk_summary.csv", index=False, float_format="%.6f")

    print(f"\nwrote {out_dir}/fig2_ranking_per_query.csv.gz ({len(ranking)} rows)")
    print(f"wrote {out_dir}/fig2_topk_per_query.csv ({len(per_query)} queries)")
    print(f"wrote {out_dir}/fig2_topk_summary.csv ({len(summary)} rows)\n")

    above = summary[summary.above_null]
    below = summary[summary.below_null]
    print(f"methods x k with a paired 95% interval ABOVE their own null: {len(above)} "
          f"of {len(summary)}")
    print(f"methods x k with a paired 95% interval BELOW their own null: {len(below)} "
          f"of {len(summary)}")
    if len(above):
        print(above[["method", "k", "hit_rate", "null_rate", "paired_lo", "paired_hi"]]
              .to_string(index=False))


if __name__ == "__main__":
    main()
