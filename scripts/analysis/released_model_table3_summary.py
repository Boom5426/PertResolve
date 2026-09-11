#!/usr/bin/env python3
"""Reproduce the released-model entries used in Supplementary Table 3.

PDS comes from the canonical multi-seed, distinct-variant summary. Pearson-delta
is recomputed from the committed per-evaluation table by first averaging repeated
settings within each (gene, variant) pair and then averaging distinct variants.

The committed per-evaluation table records zero placeholders for scVIDR's
non-finite JAK1 predictions. They are not correlations and are excluded here.
The audit established that all and only the 59 scVIDR JAK1 evaluation rows,
covering 26 variants, are non-finite in the prediction archive. Assertions below
make that exception explicit rather than silently treating the placeholders as
measurements.

Usage:
    python scripts/analysis/released_model_table3_summary.py
    python scripts/analysis/released_model_table3_summary.py --check
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[2]
SCORES = REPO / "results" / "canonical" / "unified_results5.csv"
PDS = REPO / "results" / "canonical" / "definitive_summary.csv"
OUTPUT = REPO / "results" / "canonical" / "released_model_table3_summary.csv"
METHODS = ("CellFlow", "Biolord", "scGen", "PerturbNet", "scVIDR")


def compute() -> pd.DataFrame:
    scores = pd.read_csv(SCORES)
    pds = pd.read_csv(PDS).set_index("method")
    rows: list[dict[str, object]] = []

    for method in METHODS:
        block = scores.loc[scores["method"].eq(method)].copy()
        if method == "scVIDR":
            refused = block.loc[block["gene"].eq("JAK1")]
            assert len(refused) == 59
            assert refused[["gene", "variant"]].drop_duplicates().shape[0] == 26
            assert np.allclose(refused["pearson_delta"], 0.0)
            block = block.loc[block["gene"].ne("JAK1")]

        finite = block.loc[np.isfinite(block["pearson_delta"])]
        per_variant = (
            finite.groupby(["gene", "variant"], as_index=False)["pearson_delta"]
            .mean()
        )
        canonical = pds.loc[method]
        rows.append(
            {
                "method": method,
                "pds_query_n": int(canonical["n_variants"]),
                "pds": float(canonical["PDS"]),
                "pds_ci_lo": float(canonical["ci_lo"]),
                "pds_ci_hi": float(canonical["ci_hi"]),
                "pds_refused_variants": int(canonical["n_refused"]),
                "pearson_finite_rows": len(finite),
                "pearson_distinct_variants": len(per_variant),
                "pearson_delta": float(per_variant["pearson_delta"].mean()),
            }
        )

    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="compare the recomputed table with the committed output",
    )
    args = parser.parse_args()
    result = compute()

    if args.check:
        expected = pd.read_csv(OUTPUT)
        pd.testing.assert_frame_equal(
            result,
            expected,
            check_exact=False,
            rtol=1e-12,
            atol=1e-12,
        )
        print(f"PASS {OUTPUT.relative_to(REPO)}")
        return

    result.to_csv(OUTPUT, index=False, float_format="%.12g")
    print(f"wrote {OUTPUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
