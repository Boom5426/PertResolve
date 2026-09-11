#!/usr/bin/env python3
"""Rebuild ``summary.csv``, the nine points Fig. 5c plots, from the committed JSONs.

The table had no generator. README.md said so, which meant the one file standing between
the per-dataset runs and a main-text panel was hand-assembled and unverifiable: nothing
could tell you whether a number in it still matched the JSON it came from.

Everything in the table except ``level`` is a rename of a JSON field. ``level`` is the
perturbation class the screen actually manipulates, which is an editorial judgement about
the dataset rather than an output of the run, so it is declared here and asserted to cover
every JSON present.

Row order is preserved rather than sorted. Fig. 5c draws markers with ``iterrows`` at one
zorder, so for the overlapping allele points the file order decides which is on top; a
re-sort would silently redraw the panel.

Usage:
    python results/benchmark_resolution/make_summary.py --out OUT
    python results/benchmark_resolution/make_summary.py --out OUT --check

    --out     directory that receives summary.csv (required, no default).
    --check   additionally diff the result against the committed summary.csv and exit
              non-zero on any disagreement.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import reject_repo_results, repo_root, require_inputs  # noqa: E402

#: What each screen perturbs, and the order Fig. 5c draws them in. Declared rather than
#: derived: no field of the run records whether a perturbation is a knockdown, a gene
#: pair, a drug or an allele.
LEVELS: dict[str, str] = {
    "Adamson": "gene-KD",
    "Norman": "gene-pair",
    "Replogle": "gene-KD",
    "VCC": "gene-KD",
    "allele_JAK1": "allele",
    "sciPlex": "drug",
    "allele_GATA1": "allele",
    "allele_KRAS": "allele",
    "allele_TP53": "allele",
}

#: JSON field -> column name in the table Fig. 5c reads.
RENAMES = {
    "n_pert": "n_perturbations",
    "frac_rankable": "rankable_fraction",
    "reliability_P_recover_order": "resolution_P_recover_order",
}
COLUMNS = ["dataset", "level", "n_perturbations", "rankable_fraction",
           "oracle_ceiling", "resolution_P_recover_order"]


def build(bench_dir: Path) -> pd.DataFrame:
    """One row per committed JSON, in the declared order."""
    found = {json.loads(p.read_text())["dataset"]: json.loads(p.read_text())
             for p in sorted(bench_dir.glob("*.json"))}

    undeclared = sorted(set(found) - set(LEVELS))
    missing = sorted(set(LEVELS) - set(found))
    if undeclared:
        raise SystemExit(
            f"these runs have no declared perturbation level: {undeclared}. Add them to "
            "LEVELS rather than letting the panel plot an unlabelled point.")
    if missing:
        raise SystemExit(
            f"these datasets are declared but have no JSON in {bench_dir}: {missing}. "
            "Re-run them, or remove them from LEVELS.")

    rows = []
    for name in LEVELS:                          # declared order, see module docstring
        d = found[name]
        row = {"dataset": name, "level": LEVELS[name]}
        row.update({new: d[old] for old, new in RENAMES.items()})
        row["oracle_ceiling"] = d["oracle_ceiling"]
        rows.append(row)
    return pd.DataFrame(rows)[COLUMNS]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, type=Path,
                    help="Directory that receives summary.csv. Required and never "
                         "defaulted, so a re-run cannot overwrite the committed table.")
    ap.add_argument("--check", action="store_true",
                    help="diff the rebuilt table against the committed summary.csv")
    args = ap.parse_args()
    out_dir = reject_repo_results(args.out)

    bench_dir = repo_root() / "results" / "benchmark_resolution"
    require_inputs(bench_dir)
    df = build(bench_dir)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "summary.csv"
    df.to_csv(out_path, index=False)
    print(df.to_string(index=False))
    print(f"\nwrote {out_path}")

    if args.check:
        committed = bench_dir / "summary.csv"
        require_inputs(committed)
        ref = pd.read_csv(committed)
        if ref.shape != df.shape or list(ref.columns) != list(df.columns):
            raise SystemExit(f"shape or columns differ: committed {ref.shape} "
                             f"{list(ref.columns)} vs rebuilt {df.shape} "
                             f"{list(df.columns)}")
        same = ref.reset_index(drop=True).equals(df.reset_index(drop=True))
        if not same:
            diff = ref.compare(df, result_names=("committed", "rebuilt"))
            raise SystemExit(f"rebuilt table disagrees with the committed one:\n{diff}")
        print(f"check: rebuilt table is identical to {committed}")


if __name__ == "__main__":
    main()
