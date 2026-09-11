#!/usr/bin/env python3
"""Re-pool the released per-seed scores, refusing predictions that were never scorable.

``definitive_summary.csv`` is produced by ``train_only_grid.py``, which needs the
per-gene cell arrays. Those are not redistributed, so the summary could not be
regenerated from a fresh checkout even though every number in it is a pooling of
``results/canonical/per_seed_variant_pds.csv.gz``, which is released. This script
closes that gap: it recomputes the whole table from the released per-seed scores
and reproduces the committed values exactly, then applies one documented refusal.

The refusal exists because the harness scores an unscorable prediction at the
chance value. ``harness.py`` returns 0.5 when the distance to the true match is
non-finite, so scVIDR, whose re-implementation diverged on JAK1 and returned no
finite prediction for any of that gene's 26 variants, carries 59 held-out
evaluations that are an imputation rather than a measurement. Averaging them into
a reported estimate mixes the two, so they are refused and the coverage is
reported instead.

Refusal is by explicit (method, gene) pair, not by detection. A rule such as "all
per-seed values are exactly 0.5" would also catch the WT-null reference, whose
0.5 is the designed consequence of predicting a zero vector and is a measurement
of that reference, not a failure. The two look identical in the table and only the
provenance separates them, so the pair is named here and cross-checked against
``n_nonfinite`` in ``results/canonical/residual_axis_models.csv``.

Usage:
    python scripts/analysis/repool_per_seed.py --check
    python scripts/analysis/repool_per_seed.py --out DIR [--convention refuse]

``--check`` recomputes under the imputation convention and compares against the
committed table without writing, which is what pins this script to the published
run. ``--out`` writes ``definitive_summary.csv``; the directory may not be inside
``results/``, matching the other generators here.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import reject_repo_results, repo_root, require_inputs

NBOOT = 2000
SEED = 0
MODE = "train-only"

# (method, gene) cells whose stored score is the harness's chance value standing in
# for a prediction that could not be ranked at all. See the module docstring for
# why this is a named list rather than a detection rule.
REFUSED: tuple[tuple[str, str], ...] = (("scVIDR", "JAK1"),)

CHANCE = 0.50
# Every column of the committed table, so --check compares all of them. An earlier
# version listed only the numeric fields and silently dropped the two booleans,
# which the figure panels read by name.
COLUMNS = ["PDS", "ci_lo", "ci_hi", "crosses", "n_variants",
           "PDS_cell", "cell_lo", "cell_hi", "cell_crosses", "n_cells"]


def pool(sub: pd.DataFrame) -> dict[str, float | int]:
    """Pool one method's seed-averaged scores on both resampling units.

    Reproduces the pooling of ``train_only_grid.py``: the primary estimate averages
    distinct (gene, variant) pairs, so an allele held out by several splits counts
    once; the secondary estimate averages within each (gene, split) cell.

    Args:
        sub: seed-averaged rows for one method, with ``gene``, ``variant``,
            ``split`` and ``pds`` columns.

    Returns:
        The eight summary fields, rounded as the committed table rounds them.
    """
    rng = np.random.default_rng(SEED)
    per_var = sub.groupby(["gene", "variant"])["pds"].mean().to_numpy()
    boot_v = np.array([per_var[rng.integers(0, len(per_var), len(per_var))].mean()
                       for _ in range(NBOOT)])
    lo_v, hi_v = np.percentile(boot_v, [2.5, 97.5])

    cells = [g["pds"].to_numpy() for _, g in sub.groupby(["gene", "split"])]
    rng = np.random.default_rng(SEED)
    boot_c = np.array([np.mean([c[rng.integers(0, len(c), len(c))].mean() for c in cells])
                       for _ in range(NBOOT)])
    lo_c, hi_c = np.percentile(boot_c, [2.5, 97.5])

    return dict(
        PDS=round(float(per_var.mean()), 4),
        ci_lo=round(float(lo_v), 4),
        ci_hi=round(float(hi_v), 4),
        crosses=bool(lo_v <= CHANCE <= hi_v),
        n_variants=int(sub.groupby(["gene", "variant"]).ngroups),
        PDS_cell=round(float(np.mean([c.mean() for c in cells])), 4),
        cell_lo=round(float(lo_c), 4),
        cell_hi=round(float(hi_c), 4),
        cell_crosses=bool(lo_c <= CHANCE <= hi_c),
        n_cells=len(cells),
    )


def summarize(seed_avg: pd.DataFrame, refuse: bool) -> pd.DataFrame:
    """Pool every method, optionally dropping the refused (method, gene) cells."""
    rows = []
    for method, sub in seed_avg.groupby("method"):
        dropped = 0
        if refuse:
            mask = np.zeros(len(sub), bool)
            for refused_method, refused_gene in REFUSED:
                if method == refused_method:
                    mask |= (sub.gene == refused_gene).to_numpy()
            dropped = int(sub.loc[mask].groupby(["gene", "variant"]).ngroups)
            sub = sub.loc[~mask]
            if sub.empty:
                raise SystemExit(f"refusing every cell of {method} leaves nothing to pool")
        rows.append(dict(method=method, standardization=MODE, **pool(sub),
                         n_refused=dropped))
    return pd.DataFrame(rows).sort_values("PDS", ascending=False).reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", default=None,
                    help="per_seed_variant_pds.csv.gz "
                         "(default: results/canonical/per_seed_variant_pds.csv.gz)")
    ap.add_argument("--convention", choices=("impute", "refuse"), default="refuse",
                    help="impute reproduces the published table; refuse drops the "
                         "unscorable cells and reports coverage (default)")
    ap.add_argument("--out", default=None, metavar="OUT_DIR",
                    help="directory receiving definitive_summary.csv; may not be "
                         "inside results/")
    ap.add_argument("--check", action="store_true",
                    help="recompute under the impute convention and assert it "
                         "matches the committed table, then exit")
    args = ap.parse_args()

    scores_path = Path(args.scores) if args.scores else \
        repo_root() / "results" / "canonical" / "per_seed_variant_pds.csv.gz"
    require_inputs(scores_path)
    long = pd.read_csv(scores_path)
    seed_avg = (long.groupby(["method", "gene", "split", "variant"])["pds"]
                .mean().reset_index())

    if args.check:
        committed = (pd.read_csv(repo_root() / "results" / "canonical" /
                                 "definitive_summary.csv").set_index("method"))
        recomputed = summarize(seed_avg, refuse=False).set_index("method")
        missing = [c for c in committed.columns if c not in recomputed.columns]
        if missing:
            raise SystemExit(f"recomputed table is missing committed column(s): {missing}")
        bad = []
        for m in committed.index:
            for c in COLUMNS:
                want, got = committed.loc[m, c], recomputed.loc[m, c]
                same = (bool(want) == bool(got) if c.endswith("crosses")
                        else abs(float(want) - float(got)) <= 1e-9)
                if not same:
                    bad.append((m, c))
        if bad:
            for m, c in bad:
                print(f"MISMATCH {m}.{c}: committed "
                      f"{committed.loc[m, c]} recomputed {recomputed.loc[m, c]}")
            raise SystemExit(f"{len(bad)} field(s) differ from the committed table")
        print(f"impute convention reproduces all {len(committed)} committed rows exactly")
        return

    if not args.out:
        raise SystemExit("--out is required unless --check is given")
    out_dir = reject_repo_results(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    table = summarize(seed_avg, refuse=args.convention == "refuse")
    table.to_csv(out_dir / "definitive_summary.csv", index=False)
    print(table.to_string(index=False))
    for method, gene in REFUSED:
        row = table[table.method == method]
        if not row.empty and args.convention == "refuse":
            print(f"\n{method}: {gene} refused, {int(row.n_refused.iloc[0])} variants "
                  f"dropped, scored on {int(row.n_variants.iloc[0])}")
    print(f"\nwrote {out_dir / 'definitive_summary.csv'} "
          f"({args.convention} convention)")


if __name__ == "__main__":
    main()
