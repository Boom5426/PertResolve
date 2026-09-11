#!/usr/bin/env python3
"""Audit the design-ordered predictor family that Fig. 5c and Fig. 5d are scored on.

Two things can go wrong with a family whose quality is fixed by construction, and the
scripts that build it (``results/benchmark_resolution/benchmark_resolution.py`` and
``scripts/analysis/controlled_predictors.py``) check neither.

**Saturation.** ``pertresolve/resolution/recovery.py`` records that interpolating toward
the panel mean is degenerate when that mean is near zero: cosine PDS is scale-blind, so
every weight above zero can collapse onto one score, and the packaged path moved to a
deranged neighbour because of it. A tied ladder is indistinguishable from a wrongly
ordered one in ``P_correct_order``, because ``inversions()`` counts ``>=``. If the ladder
that Fig. 5c and Fig. 5d score is saturated, a recovery probability of 0 means "tied",
not "inverted", and any power curve built on the family measures nothing.

**Gap mislabelling.** ``gap_resolution`` in the committed JSONs is keyed on the *weight*
difference between adjacent members, not on the PDS quality difference that difference
produces. The map from weight to PDS saturates hard, so the two are not interchangeable:
a weight gap of 0.035 buys a PDS gap of about 0.002 in Replogle. Any statement of the form
"the benchmark resolves a quality gap of X" has to quote the realized PDS gap.

This script answers both from the committed tables alone; it needs no workspace.

Usage:
    python scripts/analysis/alpha_ladder_audit.py --out OUT

    --out   directory that receives alpha_ladder_audit.csv (required, no default).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import reject_repo_results, repo_root, require_inputs  # noqa: E402

#: Tie tolerance, taken from ``pertresolve.resolution.recovery.SATURATION_TOL`` so this
#: audit and the packaged guard agree on what "the same score" means.
SATURATION_TOL = 1e-3


def _ladder(pds_by_alpha: dict[str, float]) -> tuple[list[float], list[float]]:
    """Weights and their scores, sorted by weight."""
    items = sorted((float(k), float(v)) for k, v in pds_by_alpha.items())
    return [a for a, _ in items], [v for _, v in items]


def audit_external(path: Path) -> dict | None:
    """One row per external screen, or None where the ladder was not stored."""
    d = json.loads(path.read_text())
    if "pds_by_alpha" not in d:
        return None
    alphas, scores = _ladder(d["pds_by_alpha"])
    steps = list(zip(alphas, alphas[1:], scores, scores[1:]))
    ties = sum(1 for _, _, lo, hi in steps if abs(hi - lo) <= SATURATION_TOL)
    inversions = sum(1 for _, _, lo, hi in steps if hi < lo)
    tail = scores[1:]
    return dict(
        dataset=d["dataset"],
        n_pert=d["n_pert"],
        n_steps=len(steps),
        adjacent_ties=ties,
        adjacent_inversions=inversions,
        # The failure mode recovery.py names: everything above the zero weight collapses.
        saturated=bool(max(tail) - min(tail) <= SATURATION_TOL),
        pds_at_alpha0=scores[0],
        pds_at_alpha1=scores[-1],
        smallest_realized_gap=round(min(hi - lo for _, _, lo, hi in steps), 5),
        largest_realized_gap=round(max(hi - lo for _, _, lo, hi in steps), 5),
    )


def audit_steps(path: Path) -> list[dict]:
    """One row per adjacent weight pair: the weight gap, the gap it buys, and P(resolve)."""
    d = json.loads(path.read_text())
    if "pds_by_alpha" not in d:
        return []
    alphas, scores = _ladder(d["pds_by_alpha"])
    gap_resolution = d.get("gap_resolution", {})
    rows = []
    for k in range(len(alphas) - 1):
        weight_gap = round(alphas[k + 1] - alphas[k], 3)
        rows.append(dict(
            dataset=d["dataset"],
            step=k,
            alpha_lo=alphas[k],
            alpha_hi=alphas[k + 1],
            weight_gap=weight_gap,
            realized_pds_gap=round(scores[k + 1] - scores[k], 5),
            # Keyed on the weight gap in the committed JSON, which is the mislabelling.
            p_resolve=gap_resolution.get(str(weight_gap)),
        ))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, type=Path,
                    help="Directory that receives alpha_ladder_audit.csv and "
                         "alpha_ladder_steps.csv. Required and never defaulted, so an "
                         "accidental re-run cannot overwrite the committed tables under "
                         "results/.")
    out_dir = reject_repo_results(ap.parse_args().out)

    bench_dir = repo_root() / "results" / "benchmark_resolution"
    recovery_csv = repo_root() / "results" / "canonical" / "controlled_recovery.csv"
    require_inputs(bench_dir, recovery_csv)

    jsons = sorted(bench_dir.glob("*.json"))
    if not jsons:
        raise SystemExit(f"no benchmark_resolution JSONs found in {bench_dir}")

    summary = [r for r in (audit_external(p) for p in jsons) if r is not None]
    steps = [r for p in jsons for r in audit_steps(p)]
    missing = [json.loads(p.read_text())["dataset"] for p in jsons
               if "pds_by_alpha" not in json.loads(p.read_text())]

    sdf, tdf = pd.DataFrame(summary), pd.DataFrame(steps)

    print("=== ladder shape, screens whose ladder was stored ===")
    print(sdf.to_string(index=False))
    print("\n=== weight gap vs the PDS gap it actually buys ===")
    print(tdf.to_string(index=False))

    print("\n=== verdict ===")
    if sdf.saturated.any():
        print("SATURATED:", ", ".join(sdf.loc[sdf.saturated, "dataset"]))
    else:
        print("not saturated: every stored ladder separates by more than "
              f"{SATURATION_TOL} at least once above the zero weight")
    if sdf.adjacent_ties.sum() == 0 and sdf.adjacent_inversions.sum() == 0:
        print("no adjacent ties and no adjacent inversions in any stored ladder, so a "
              "P_correct_order below 1 in these datasets reflects sampling, not ties")
    else:
        print(f"adjacent ties: {int(sdf.adjacent_ties.sum())}, "
              f"adjacent inversions: {int(sdf.adjacent_inversions.sum())}")
    ratio = tdf.weight_gap / tdf.realized_pds_gap.abs()
    print(f"weight gap over realized PDS gap: median {ratio.median():.1f}x, "
          f"max {ratio.max():.1f}x -- gap_resolution keys are weight gaps, not quality gaps")

    # The four allele datasets and every row of Fig. 5d store only the ends of the ladder,
    # so the same audit cannot be run on them from committed data.
    print("\n=== not auditable from committed tables ===")
    print("benchmark_resolution JSONs without pds_by_alpha:", ", ".join(missing) or "none")
    rec = pd.read_csv(recovery_csv)
    print(f"{recovery_csv.name}: {len(rec)} rows store pds_a0 and pds_a1 only; the "
          "intermediate weights that P_correct_order is computed over are not persisted, "
          "so Fig. 5d needs a re-run to audit.")
    print("Also visible there: JAK1's cohort shrinks as depth rises "
          f"({', '.join(f'm={r.depth_m}:n={r.n_var}' for r in rec[rec.gene == 'JAK1'].itertuples())}), "
          "so depth and candidate-set size are confounded along that curve.")

    out_dir.mkdir(parents=True, exist_ok=True)
    sdf.to_csv(out_dir / "alpha_ladder_audit.csv", index=False)
    tdf.to_csv(out_dir / "alpha_ladder_steps.csv", index=False)
    print(f"\nwrote {out_dir / 'alpha_ladder_audit.csv'} and "
          f"{out_dir / 'alpha_ladder_steps.csv'}")


if __name__ == "__main__":
    main()
