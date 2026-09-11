#!/usr/bin/env python3
"""Rebuild the six-dimensional biophysical variant vector (theta) from its definition.

The committed feature table `data/pertresolve_bench.csv` has no generator: no
script in this repository or on the compute server produces it, and it entered git
in a single squashed commit. That is the root cause of two defects found in it.

  1. `is_hotspot` is not the external annotation the Methods describes. For GATA1 it
     is reproduced 255/255 by `(severity_dist >= 80th percentile) OR
     clinvar_pathogenic`, and `severity_dist` correlates with the measured
     variant-to-wild-type expression distance at Pearson 1.0000. A model input was
     therefore a thresholded readout of the response being predicted, with the
     threshold itself computed across held-out variants.
  2. The TP53 and KRAS physicochemical terms do not reproduce the amino-acid tables
     in `pertresolve/features.py`: the same substitution takes different values at
     different positions, and `d_vol` is compressed onto a near-constant pedestal.

This script rebuilds theta from `pertresolve.features.compute_theta` plus an
explicit annotation spec, so the table becomes reproducible and both defects become
detectable. It does not attempt to invert the transform that produced the committed
values; it recomputes from the definition.

It regenerates the feature columns only. The variant list, cell counts, split roles
and measured columns are passed through from `--rows` unchanged, because they are
measured or assigned quantities that no annotation spec can produce.

Usage:
  python build_bench_features.py --out /path/to/pertresolve_bench_rebuilt.csv
  python build_bench_features.py --out OUT --compare data/pertresolve_bench.csv \\
                                 --compare data/pertresolve_bench_v2.csv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pertresolve.features import compute_theta  # noqa: E402
from pertresolve.paths import (reject_repo_results, repo_root,  # noqa: E402
                                 require_inputs)

#: Column order of the feature vector, matching `harness.THETA`.
THETA_COLUMNS = ["d_hydro", "d_vol", "d_charge", "fold_core", "cat_switch", "is_hotspot"]

#: A single amino-acid substitution, 1-based position, as written in variant names.
SUBSTITUTION = re.compile(r"^([A-Z])(\d+)([A-Z])$")

#: A condition that carries a residue position but no mutant residue, for example
#: `Q387*` (stop-gained) or `syn_L158` (synonymous, named rather than spelled out).
POSITION_ONLY = re.compile(r"^(?:syn_)?[A-Z](\d+)\*?$")

#: `compute_theta` accepts `protein_length` but does not read it. Passing a sentinel
#: keeps unverified per-gene lengths out of the config.
PROTEIN_LENGTH_UNUSED = 0


class SpecError(ValueError):
    """The annotation spec is malformed or incomplete."""


class LeakageError(RuntimeError):
    """A feature column carries information about the measured response."""


def load_spec(path: Path) -> dict:
    """Load and validate the annotation spec.

    Raises:
        SpecError: if a required key is missing, a range is malformed, or the
            declared f_v tiers are not strictly decreasing in value.
    """
    with open(path) as handle:
        spec = yaml.safe_load(handle)

    for key in ("fold_core_tiers", "genes", "non_missense", "leakage_gate"):
        if key not in spec:
            raise SpecError(f"{path}: missing top-level key {key!r}")

    tiers = spec["fold_core_tiers"]
    values = [float(t["value"]) for t in tiers]
    if values != sorted(values, reverse=True):
        raise SpecError(f"{path}: fold_core_tiers must be ordered most specific first, "
                        f"got values {values}")

    for gene, entry in spec["genes"].items():
        for key in ("structural_core", "interaction_regions", "functional_switch_residues",
                    "hotspot_codons", "hotspot_regions", "hotspot_variants"):
            if key not in entry:
                raise SpecError(f"{path}: gene {gene} is missing key {key!r}")
        for key in ("structural_core", "interaction_regions", "hotspot_regions"):
            for rng in entry[key]:
                if len(rng) != 2 or int(rng[0]) > int(rng[1]):
                    raise SpecError(f"{path}: gene {gene}, {key}: malformed range {rng}")
    return spec


def parse_variant(name: str, wildtype_labels: tuple[str, ...]) -> tuple[str, list]:
    """Classify a condition name and return its substitutions or positions.

    Returns:
        ``(kind, payload)`` where ``kind`` is one of ``"wildtype"``, ``"missense"``
        or ``"position_only"``. For ``"missense"`` the payload is a list of
        ``(wt_aa, position, mut_aa)``; for ``"position_only"`` a list of positions,
        possibly empty when the name carries none, as for a splice condition.
    """
    text = str(name).strip()
    if text in wildtype_labels:
        return "wildtype", []

    parts = [p.strip() for p in text.split(",")]
    matches = [SUBSTITUTION.match(p) for p in parts]
    if all(matches):
        return "missense", [(m.group(1), int(m.group(2)), m.group(3)) for m in matches]

    positions = []
    for part in parts:
        hit = POSITION_ONLY.match(part)
        if hit:
            positions.append(int(hit.group(1)))
    return "position_only", positions


def in_any(position: int, ranges: list) -> bool:
    """True when a 1-based position falls inside any inclusive range."""
    return any(int(lo) <= position <= int(hi) for lo, hi in ranges)


def fold_core_value(positions: list[int], entry: dict, tiers: list[dict]) -> float:
    """Return f_v for a condition, taking the first tier any mutated position hits."""
    by_name = {t["name"]: float(t["value"]) for t in tiers}
    for tier in tiers:
        key = {"structural_core": "structural_core",
               "interaction_region": "interaction_regions"}[tier["name"]]
        if any(in_any(p, entry[key]) for p in positions):
            return by_name[tier["name"]]
    return 0.0


def annotations(name: str, positions: list[int], entry: dict, tiers: list[dict]) -> tuple:
    """Return ``(f_v, s_v, h_v)`` for one condition from the spec."""
    f_v = fold_core_value(positions, entry, tiers)
    switch = set(int(r) for r in entry["functional_switch_residues"])
    s_v = 1.0 if any(p in switch for p in positions) else 0.0
    codons = set(int(c) for c in entry["hotspot_codons"])
    h_v = 1.0 if (any(p in codons for p in positions)
                  or any(in_any(p, entry["hotspot_regions"]) for p in positions)
                  or str(name).strip() in set(entry["hotspot_variants"])) else 0.0
    return f_v, s_v, h_v


def build_row(name: str, gene: str, spec: dict) -> np.ndarray:
    """Compute the six-vector for one condition."""
    entry = spec["genes"][gene]
    tiers = spec["fold_core_tiers"]
    wt_labels = tuple(spec["non_missense"]["wildtype_labels"])
    kind, payload = parse_variant(name, wt_labels)

    if kind == "wildtype":
        return np.zeros(6, dtype=float)

    if kind == "missense":
        positions = [p for _, p, _ in payload]
        # compute_theta is authoritative for the physicochemical terms; annotations
        # are recomputed here so the declared f_v tiers can be applied.
        per_sub = np.array([
            compute_theta(wt, mut, pos, PROTEIN_LENGTH_UNUSED)[:3]
            for wt, pos, mut in payload
        ], dtype=float)
        physchem = per_sub.mean(axis=0)
    else:
        positions = payload
        if spec["non_missense"]["physicochemical"] != "zero":
            raise SpecError("non_missense.physicochemical: only 'zero' is implemented")
        physchem = np.zeros(3, dtype=float)
        if spec["non_missense"]["annotations"] != "by_position":
            raise SpecError("non_missense.annotations: only 'by_position' is implemented")

    return np.concatenate([physchem, np.array(annotations(name, positions, entry, tiers))])


def cross_check_binary_fold_core(frame: pd.DataFrame, spec: dict) -> None:
    """Assert the rebuilt f_v equals `compute_theta`'s when no 0.5 tier is declared.

    Keeps `pertresolve.features` the single source of truth for the binary case,
    so a divergence means this script and the package have drifted apart.
    """
    wt_labels = tuple(spec["non_missense"]["wildtype_labels"])
    for gene, entry in spec["genes"].items():
        if entry["interaction_regions"]:
            continue
        sub = frame[frame.gene == gene]
        for _, row in sub.iterrows():
            kind, payload = parse_variant(row["variant"], wt_labels)
            if kind != "missense":
                continue
            ranges = [tuple(map(int, r)) for r in entry["structural_core"]]
            expected = max(
                compute_theta(wt, mut, pos, PROTEIN_LENGTH_UNUSED, domain_ranges=ranges)[3]
                for wt, pos, mut in payload
            )
            if not np.isclose(expected, row["fold_core"]):
                raise SpecError(
                    f"fold_core disagrees with pertresolve.features for "
                    f"{gene} {row['variant']}: package {expected}, rebuilt {row['fold_core']}")


def check_tiers(frame: pd.DataFrame, spec: dict) -> None:
    """Assert every annotation column only takes values the spec declares legal."""
    legal_f = {0.0} | {float(t["value"]) for t in spec["fold_core_tiers"]}
    for column, legal in (("fold_core", legal_f), ("cat_switch", {0.0, 1.0}),
                          ("is_hotspot", {0.0, 1.0})):
        seen = set(np.round(frame[column].astype(float), 9))
        illegal = seen - legal
        if illegal:
            raise SpecError(f"{column} takes undeclared values {sorted(illegal)}; "
                            f"declared {sorted(legal)}")


def response_correlation(frame: pd.DataFrame, spec: dict) -> pd.DataFrame:
    """Report how strongly each feature tracks the measured per-variant response.

    This is a diagnostic, not a gate. A feature computed without ever reading the
    response can still correlate with it, and a useful one should: variants in a
    structural core really do have larger effects. The rebuilt table reaches
    Spearman +0.41 for KRAS `d_vol` with no leakage whatsoever. Failing on
    correlation alone would therefore reject correct features.

    The property that actually distinguishes a leak is provenance, and provenance is
    enforced by :func:`audit`: a column that cannot be regenerated from the variant
    name and the spec is carrying information from somewhere else. That is the
    invariant `is_hotspot` violated.
    """
    gate = spec["leakage_gate"]
    column = gate["leak_column"]
    if column not in frame.columns:
        raise LeakageError(f"leak column {column!r} is not in the row source, so the "
                           f"diagnostic cannot run; supply it or change leakage_gate")

    def spearman(a: np.ndarray, b: np.ndarray) -> float:
        return float(np.corrcoef(np.argsort(np.argsort(a)), np.argsort(np.argsort(b)))[0, 1])

    rows = []
    for gene in sorted(frame.gene.unique()):
        sub = frame[(frame.gene == gene) & frame[column].notna()]
        if len(sub) < int(gate["min_variants"]):
            continue
        for feature in THETA_COLUMNS:
            values = sub[feature].astype(float).values
            rho = 0.0 if np.std(values) == 0 else spearman(
                values, sub[column].astype(float).values)
            rows.append((gene, feature, len(sub), rho))
    return pd.DataFrame(rows, columns=["gene", "feature", "n", "spearman_vs_measured"])


def audit(candidate: Path, spec: dict, rows_source: pd.DataFrame) -> pd.DataFrame:
    """Fail if a candidate table's feature columns are not reproducible from the spec.

    This is the hard gate. Run against `data/pertresolve_bench.csv` it fails on
    the two known defects; run against a table this script produced it passes.

    Raises:
        LeakageError: if any feature column of the candidate differs from the rebuild.
    """
    ref = pd.read_csv(candidate)
    built = np.array([build_row(row["variant"], row["gene"], spec)
                      for _, row in rows_source.iterrows()])
    rebuilt = rows_source[["gene", "variant"]].copy()
    for index, column in enumerate(THETA_COLUMNS):
        rebuilt[column] = built[:, index]

    merged = rebuilt.merge(ref[["gene", "variant"] + THETA_COLUMNS],
                           on=["gene", "variant"], suffixes=("_new", "_ref"))
    rows, failures = [], []
    for gene in sorted(merged.gene.unique()):
        sub = merged[merged.gene == gene]
        for column in THETA_COLUMNS:
            differing = ~np.isclose(sub[f"{column}_new"].astype(float),
                                    sub[f"{column}_ref"].astype(float), atol=1e-9)
            rows.append((gene, column, int(differing.sum()), len(sub)))
            if differing.any():
                failures.append(f"{gene}.{column}: {int(differing.sum())} of {len(sub)} "
                                f"rows are not reproducible from the spec")
    report = pd.DataFrame(rows, columns=["gene", "column", "n_differ", "n_rows"])
    if failures:
        raise LeakageError(
            f"{candidate} is not reproducible from {spec.get('_path', 'the spec')}:\n  "
            + "\n  ".join(failures))
    return report


def compare(frame: pd.DataFrame, reference: Path) -> pd.DataFrame:
    """Per gene and column, count rows where the rebuild differs from a reference."""
    ref = pd.read_csv(reference)
    merged = frame[["gene", "variant"] + THETA_COLUMNS].merge(
        ref[["gene", "variant"] + THETA_COLUMNS], on=["gene", "variant"],
        suffixes=("_new", "_ref"))
    rows = []
    for gene in sorted(merged.gene.unique()):
        sub = merged[merged.gene == gene]
        for column in THETA_COLUMNS:
            differing = ~np.isclose(sub[f"{column}_new"].astype(float),
                                    sub[f"{column}_ref"].astype(float), atol=1e-9)
            rows.append((gene, column, int(differing.sum()), len(sub)))
    return pd.DataFrame(rows, columns=["gene", "column", "n_differ", "n_rows"])


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rows", type=Path,
                        default=repo_root() / "data" / "pertresolve_bench.csv",
                        help="table supplying the condition list and the passthrough "
                             "columns (default: %(default)s)")
    parser.add_argument("--spec", type=Path,
                        default=repo_root() / "configs" / "bench_features.yaml",
                        help="annotation spec (default: %(default)s)")
    parser.add_argument("--out", type=Path,
                        help="destination CSV for the rebuilt table; omit with --audit")
    parser.add_argument("--compare", type=Path, action="append", default=[],
                        help="existing table to diff the rebuild against; repeatable")
    parser.add_argument("--audit", type=Path, action="append", default=[],
                        help="fail unless this table's feature columns are exactly "
                             "reproducible from the spec; repeatable. This is the gate")
    parser.add_argument("--force", action="store_true",
                        help="permit writing over a table tracked in data/")
    args = parser.parse_args()
    if args.out is None and not args.audit:
        parser.error("give --out to build, or --audit to check an existing table")

    require_inputs(args.rows, args.spec)
    if args.out is not None:
        reject_repo_results(args.out)
    tracked = (repo_root() / "data").resolve()
    if args.out is not None and args.out.resolve().parent == tracked and not args.force:
        raise SystemExit(f"refusing to write into {tracked} without --force; the "
                         f"canonical tables are inputs to every published number")

    spec = load_spec(args.spec)
    spec["_path"] = str(args.spec)
    source = pd.read_csv(args.rows)
    missing = set(source.gene.unique()) - set(spec["genes"])
    if missing:
        raise SpecError(f"{args.spec} has no entry for gene(s) {sorted(missing)}")

    built = np.array([build_row(row["variant"], row["gene"], spec)
                      for _, row in source.iterrows()])
    out = source.copy()
    for index, column in enumerate(THETA_COLUMNS):
        out[column] = built[:, index]
    for column in ("d_hydro_raw", "d_vol_raw", "d_charge_raw"):
        if column in out.columns:
            out[column] = built[:, THETA_COLUMNS.index(column.replace("_raw", ""))]

    check_tiers(out, spec)
    cross_check_binary_fold_core(out, spec)

    print(f"rebuilt {len(out)} conditions from {args.spec}")
    corr = response_correlation(out, spec)
    limit = float(spec["leakage_gate"]["max_abs_spearman"])
    flagged = corr[corr.spearman_vs_measured.abs() > limit]
    print("\n=== diagnostic: Spearman of each rebuilt feature against the measured "
          "response ===")
    print("  (advisory only; a correct feature may legitimately correlate with the "
          "outcome)")
    print(flagged.to_string(index=False) if len(flagged)
          else f"  nothing above |{limit}|")

    for reference in args.compare:
        require_inputs(reference)
        print(f"\n=== rows differing from {reference} ===")
        table = compare(out, reference)
        print(table[table.n_differ > 0].to_string(index=False)
              if (table.n_differ > 0).any() else "  identical on every column")

    for candidate in args.audit:
        require_inputs(candidate)
        print(f"\n=== audit: is {candidate} reproducible from the spec? ===")
        audit(candidate, spec, source)
        print("  yes, every feature column matches")

    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.out, index=False)
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
