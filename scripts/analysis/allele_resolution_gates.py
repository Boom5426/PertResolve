#!/usr/bin/env python3
"""Detection and local-competition measurements for the four allele genes, from one measurement.

The historical Figure 4 analysis read its detection quantity and its nearest-competitor quantity from two
different pipelines. ``results/second_probe_rankability_table.csv`` measures detection with
an energy distance over 50 random halvings of a cell pool truncated at 300 cells per variant;
``results/canonical/resolution_sweep.csv`` measures nearest-competitor separation with a
debiased squared Euclidean distance over 10 seeds of four disjoint groups cut from the full
pool. The two can be aligned on gene, on nominal depth and on the variant list, but never on
the cell pool, the split construction or the distance, so a panel plotting one against the
other would be a join across three unmatched axes. This script removes the seam instead of
documenting it: :func:`pertresolve.resolution.resolution_report` cuts the cells once and
derives both gates from that single cut, so every point it produces is internally matched.

The two axes it reports are the released criteria, not proxies for them. They remain
independent measurements; detection is not treated as a logical prerequisite for
identification.

**Detection.** ``signal > width``, where ``signal = d_null - d_self`` is how far a variant
sits from wild type beyond its own split-half reproducibility floor, and ``width`` is the 95%
spread of ``d_self`` across splits (``pertresolve/resolution/window.py``). The continuous
form reported here is ``detect_margin = signal / width``, which crosses its gate at 1.

**Local competition.** ``separation > spread``, where ``separation`` is the cross-fitted,
debiased nearest-competitor squared separation averaged over splits and ``spread`` is its own
95% spread across those same splits (``pertresolve/resolution/report.py``). The continuous
form is ``ident_margin = separation / spread``, which also crosses its gate at 1.

Both margins are gated at 1 because both criteria compare a signal with the 95% spread of its
own sampling distribution. That is the point of running them together: the two axes of the
resulting panel are the same kind of quantity, cut from one shuffle of one cell pool.

They are not, however, measured on the same cells, and the distinction is worth stating
because it is the stronger property. At a given seed :func:`group_profiles` cuts four disjoint
blocks per variant, and :func:`detection_window` reproduces blocks 0 and 1 exactly, while the
nearest-competitor separations are formed from blocks 2 and 3. So the detection estimate and
the competition estimate share a cell pool, a control reference, a depth and a seed stream,
and share no cell. Neither gate can inherit the other's noise realisation.

**``R = d_self / d_null`` is reported but is not a gate.** On the committed rankability table
963 of 1,664 allele rows have ``R < 1`` while only 50 clear ``signal > width``, so a panel
that drew ``R = 1`` as a detection boundary would overstate the detectable share by roughly
nineteen-fold. ``R`` is kept in the per-unit output because Fig. 4a defines it and Fig. 4c
plots it, and because the distance between ``R = 1`` and the real gate is itself worth
reading. It must not be relabelled as the criterion.

Two representations are run, and neither is privileged.

``pca``  50 components, which is the frozen resolution-panel criterion (depth 50, n_seeds 8,
         n_boot 1000, seed 0, n_components 50) that all 31 external configurations in Figure 5
         were measured on. Running the allele genes on it makes them directly comparable to
         those screens for the first time.
``raw``  the gene matrices as ``harness.load_gene`` returns them, which is the space
         ``resolution_sweep.py`` measured in and therefore the space Fig. 4e reads.

If the two disagree on an axis result, that disagreement is the result and belongs in the
panel's caption; it is not a reason to pick the friendlier one.

**What may and may not be joined to results/canonical/resolution_sweep.csv.** ``eta2``,
``rho2`` and ``rho2_nn_median`` are invariant to the wild-type reference, which cancels in
both terms, so they differ from the sweep only through the seed stream and the split count and
are broadly comparable. ``split_half_reproducibility_reference``, ``p_correct_order`` and
``p_correct_winner``
are not: the sweep gives its query and truth groups two disjoint wild-type half-means while
:func:`group_profiles` shares one reference across all four blocks, and the two use different
seed streams and different split counts. Joining those three to the sweep's ``ceiling_pds`` or
``p_correct`` on (gene, depth) would silently compare two estimators. The columns are kept
because they came free with the report, not because they are cross-comparable.

Cohort. Variants are restricted to the committed benchmark table's non-wild-type rows for the
gene, which is ``resolution_sweep.py``'s rule, and :func:`group_profiles` then applies the
``4 * depth`` cell requirement itself. The resulting counts reproduce the ``real`` rows of
``resolution_sweep.csv`` (98/98/68, 92/92/84, 247/209/107, 14/9/5) and the script asserts
this, so a cohort that has silently moved fails the run rather than the reader.

The PCA is fitted once per gene on every cell of that gene, wild type included, which is what
``scripts/analysis/resolution_panel_v2.py`` does for the external screens. It is unsupervised
and sees no split, but it does see all the cells a later split will be drawn from; the choice
is inherited from the frozen criterion rather than made here, and changing it would break
comparability with Figure 5.

Usage:
  python allele_resolution_gates.py --out OUT_DIR [--base BASE] [--spaces pca,raw]
                                    [--depths 25,50,100] [--genes TP53,KRAS,GATA1,JAK1]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import (
    add_harness_to_path,
    reject_repo_results,
    require_inputs,
    resolve_base,
)
from pertresolve.resolution import resolution_report

#: The benchmark table the committed model grid was fitted on. Two tables exist and they
#: differ on 182 of 472 rows (docs/POST_FREEZE_ISSUES.md, P1-01); picking the wrong one
#: silently changes the cohort, so the digest is checked rather than the path trusted.
CANONICAL_BENCH_MD5 = "0c4869f955f78e0ab0c6e71df1ad3274"

#: Cohort sizes of the ``real`` rows of results/canonical/resolution_sweep.csv, which this
#: run must reproduce for its points to be comparable with Fig. 4e. Keyed (gene, depth).
EXPECTED_N: dict[tuple[str, int], int] = {
    ("TP53", 25): 98, ("TP53", 50): 98, ("TP53", 100): 68,
    ("KRAS", 25): 92, ("KRAS", 50): 92, ("KRAS", 100): 84,
    ("GATA1", 25): 247, ("GATA1", 50): 209, ("GATA1", 100): 107,
    ("JAK1", 25): 14, ("JAK1", 50): 9, ("JAK1", 100): 5,
}

WT_TAGS = ("WT", "wt", "WT_control")

_ap = argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--base", default=None,
                 help="processed-data directory holding the gene arrays "
                      "(env: PERTRESOLVE_DATA)")
_ap.add_argument("--out", required=True, metavar="OUT_DIR",
                 help="directory receiving the two tables; may not be inside the "
                      "repository's results/")
_ap.add_argument("--genes", default="TP53,KRAS,GATA1,JAK1",
                 help="comma-separated genes (default: %(default)s)")
_ap.add_argument("--depths", default="25,50,100",
                 help="comma-separated cells per group (default: %(default)s)")
_ap.add_argument("--spaces", default="pca,raw",
                 help="representations to run: pca, raw, or both (default: %(default)s)")
_ap.add_argument("--n-seeds", type=int, default=8,
                 help="cell splits; the frozen panel criterion is 8 (default: %(default)s)")
_ap.add_argument("--seed", type=int, default=0, help="base seed (default: %(default)s)")
_ap.add_argument("--n-boot", type=int, default=1000,
                 help="bootstraps for ordering recovery (default: %(default)s)")
_ap.add_argument("--n-components", type=int, default=50,
                 help="PCA components for the pca space (default: %(default)s)")
_ap.add_argument("--allow-cohort-drift", action="store_true",
                 help="continue when a cohort size no longer matches resolution_sweep.csv. "
                      "Off by default: a moved cohort means the new points are not "
                      "comparable with Fig. 4e and the caller must say so deliberately.")
_args = _ap.parse_args()

BASE = resolve_base(_args.base)
OUT_DIR = reject_repo_results(_args.out)
OUT_DIR.mkdir(parents=True, exist_ok=True)

BENCH_PATH = BASE / "pertresolve_bench.csv"
require_inputs(BENCH_PATH)
_digest = hashlib.md5(BENCH_PATH.read_bytes()).hexdigest()
if _digest != CANONICAL_BENCH_MD5:
    raise SystemExit(
        f"{BENCH_PATH} has md5 {_digest}, not the canonical {CANONICAL_BENCH_MD5}. "
        "Two benchmark tables exist and they differ on 182 of 472 rows; this run would "
        "measure a different cohort from every committed table. See "
        "docs/POST_FREEZE_ISSUES.md, P1-01.")

add_harness_to_path(BASE)
import harness as H  # noqa: E402

GENES = [g.strip() for g in _args.genes.split(",") if g.strip()]
DEPTHS = [int(x) for x in _args.depths.split(",")]
SPACES = [s.strip() for s in _args.spaces.split(",") if s.strip()]
for _s in SPACES:
    if _s not in ("pca", "raw"):
        raise SystemExit(f"unknown space {_s!r}; choose from pca, raw")

_BENCH = pd.read_csv(BENCH_PATH)
BENCH = {g: set(_BENCH[_BENCH.gene == g]["variant"]) - set(WT_TAGS) for g in GENES}


def control_label(labels: np.ndarray) -> str:
    """The single wild-type tag present in this gene's labels.

    Raises rather than guessing: ``resolution_report`` takes one control string, and a gene
    carrying two wild-type spellings would have half its control cells silently scored as a
    perturbation.
    """
    present = [t for t in WT_TAGS if int((labels == t).sum()) > 0]
    if len(present) != 1:
        raise SystemExit(f"expected exactly one wild-type tag, found {present}")
    return present[0]


def reduce_pca(X: np.ndarray, n_components: int, seed: int) -> np.ndarray:
    """Dense PCA, the ``pca`` path of scripts/analysis/resolution_panel_v2.py::_reduce."""
    from sklearn.decomposition import PCA
    if n_components >= X.shape[1]:
        raise SystemExit(f"n_components {n_components} is not below the {X.shape[1]} "
                         "features; the reduction would be a no-op under a name")
    return PCA(n_components=n_components, random_state=seed).fit_transform(
        np.asarray(X, dtype=np.float64))


def margin(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    """``numerator / denominator``, NaN where the denominator is not strictly positive.

    Both gates compare a signal with the 95% spread of its own sampling distribution. A
    non-positive spread means the spread was not estimable, and a ratio formed against it
    would print as a number while carrying no information, so it is withheld. The boolean
    boolean axis result for such a unit still comes from the library's own comparison, never from this
    ratio, so a withheld margin never silently changes a gate.
    """
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    out = np.full(numerator.shape, np.nan)
    ok = denominator > 0
    out[ok] = numerator[ok] / denominator[ok]
    return out


def fraction_below_one(values: np.ndarray) -> float:
    """Share of estimable values below 1, with unestimable ones excluded from the denominator.

    ``np.asarray([nan, 0.5, 1.2]) < 1.0`` is ``[False, True, False]``, so a plain mean would
    count an unestimable ratio as evidence against the threshold and would divide by a
    different denominator than the median beside it in the same summary row.
    """
    values = np.asarray(values, dtype=float)
    ok = np.isfinite(values)
    return float((values[ok] < 1.0).mean()) if ok.any() else float("nan")


unit_rows: list[pd.DataFrame] = []
summary_rows: list[dict] = []
started = time.time()


def checkpoint(*, complete: bool) -> None:
    """Write both tables and the provenance file, after every configuration.

    The provenance is written on every checkpoint rather than once at the end, carrying a
    ``complete`` flag. A run that dies in the GATA1 raw configurations used to leave two CSVs
    on disk that were indistinguishable from a finished run and had no provenance beside them
    at all, which in this repository is how a partial table gets promoted by mistake.
    """
    pd.concat(unit_rows, ignore_index=True).to_csv(
        OUT_DIR / "allele_resolution_gates_units.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(
        OUT_DIR / "allele_resolution_gates_summary.csv", index=False)
    (OUT_DIR / "allele_resolution_gates_provenance.json").write_text(json.dumps(dict(
        complete=complete,
        configurations_done=[f"{r['gene']}/{r['space']}/{r['depth']}" for r in summary_rows],
        base=str(BASE), bench_md5=CANONICAL_BENCH_MD5,
        genes=GENES, depths=DEPTHS, spaces=SPACES,
        n_seeds=_args.n_seeds, seed=_args.seed, n_boot=_args.n_boot,
        n_components=_args.n_components,
        detection_gate="signal > width, signal = d_null - d_self, "
                       "width = 95% spread of d_self across splits",
        competition_gate="separation > spread, spread = 95% spread of the nearest-competitor "
                         "separation across the same splits",
        note="R = d_self / d_null is reported but is not a gate. "
             "split_half_reproducibility_reference, "
             "p_correct_order and p_correct_winner must not be joined to "
             "results/canonical/resolution_sweep.csv; see the module docstring.",
        numpy=np.__version__, pandas=pd.__version__,
        seconds_total=round(time.time() - started, 1),
    ), indent=2) + "\n")

for gene in GENES:
    X_raw, labels = H.load_gene(gene)
    labels = np.asarray(labels)
    control = control_label(labels)
    n_cells_by_variant = pd.Series(labels).value_counts().to_dict()
    wanted = sorted(BENCH[gene])

    # One representation is materialised at a time and released before the next. GATA1 is
    # 149,161 x 2,477, so a float64 copy is 3.0 GB; resolution_profiles makes one of its own
    # from whatever it is handed, and holding a second here doubles that for no benefit. The
    # raw matrix is therefore passed at its native float32 and upcast once, inside the
    # library, rather than twice.
    for space in SPACES:
        X = X_raw if space == "raw" else reduce_pca(X_raw, _args.n_components, _args.seed)
        for depth in DEPTHS:
            # The cohort is a pure function of the labels and 4 * depth, so it is checked
            # before the quadratic energy distances are paid for rather than after. This also
            # turns the two ways a configuration can be unmeasurable into a named skip: at
            # fewer than two eligible variants group_profiles raises, and the uncaught
            # ValueError would take the whole grid down after the earlier depths were paid.
            eligible = [v for v in wanted if int((labels == v).sum()) >= 4 * depth]
            expected = EXPECTED_N.get((gene, depth))
            if expected is not None and len(eligible) != expected:
                message = (f"{gene} depth {depth}: {len(eligible)} eligible variants, but the "
                           f"real rows of resolution_sweep.csv carry {expected}. The cohort "
                           "has moved, so these points are no longer comparable with Fig. 4e.")
                if not _args.allow_cohort_drift:
                    raise SystemExit(message + " Pass --allow-cohort-drift to proceed.")
                print("WARNING: " + message, flush=True)
            if len(eligible) < 2:
                print(f"[{gene:6} {space:3} depth {depth:3}] skipped: {len(eligible)} "
                      f"variants reach {4 * depth} cells; nothing to compare", flush=True)
                continue

            t0 = time.time()
            # This paper-specific run requests the optional cross-seed diagnostics explicitly.
            # They are recorded alongside, not substituted for, the primary axes.
            report = resolution_report(
                X, labels, control=control, depth=depth,
                n_seeds=_args.n_seeds, seed=_args.seed, n_boot=_args.n_boot,
                with_window=True, perturbations=wanted, cross_seed=True)

            window = report.window
            if window is None:
                raise SystemExit(f"{gene} {space} depth {depth}: no detection window was "
                                 "computed; the panel needs both gates from one run")
            names = list(report.details["perturbations"])
            if list(window.perturbations) != names:
                raise SystemExit(
                    f"{gene} {space} depth {depth}: the detection window covers "
                    f"{len(window.perturbations)} units and the competition side {len(names)}, "
                    "in a different order. The two gates would be paired across different "
                    "perturbations, which is the exact defect this script exists to avoid.")

            n = len(names)
            if n != len(eligible):
                raise SystemExit(f"{gene} {space} depth {depth}: the report kept {n} units "
                                 f"but the pre-check found {len(eligible)} eligible; the "
                                 "eligibility rule assumed here is not the one applied.")

            separation = np.asarray(report.details["separations"], dtype=float)
            spread = np.asarray(report.details["separation_spread"], dtype=float)
            # report.py returns identified = all-False, a placeholder rather than a measured
            # when fewer than three splits make the spread unestimable, and flags that case by
            # setting identification_fraction to NaN. Copying the placeholder through would turn
            # "cannot be judged" into a measured zero, and the jointly evaluable column would
            # inherit it as a hard 0.000 beside an honest empty cell.
            estimable = bool(np.isfinite(report.identification_fraction))
            identified = (np.asarray(report.details["identified"], dtype=bool) if estimable
                          else np.full(n, np.nan, dtype=object))

            frame = pd.DataFrame(dict(
                gene=gene, space=space, depth=depth,
                n_components=(_args.n_components if space == "pca" else X.shape[1]),
                perturbation=names,
                n_cells=[int(n_cells_by_variant.get(v, 0)) for v in names],
                # Detection, from window.py: the released gate is signal > width.
                d_self=window.d_self, d_null=window.d_null,
                signal=window.signal, width=window.width,
                detect_margin=margin(window.signal, window.width),
                detectable=window.rankable,
                # Reported for continuity with Fig. 4a/4c only. Not a gate; see the module
                # docstring for the nineteen-fold disagreement with the real criterion.
                ratio_R=window.ratio,
                R_below_one=window.ratio < 1.0,
                # Local competition, from report.py: the released gate is separation > spread.
                separation=separation, spread=spread,
                ident_margin=margin(separation, spread),
                identifiable=identified,
            ))
            frame["jointly_evaluable"] = (frame["detectable"] & frame["identifiable"]
                                          if estimable else np.nan)
            unit_rows.append(frame)

            summary = dict(gene=gene, space=space, depth=depth, n_units=n)
            summary.update(report.to_dict())
            summary.update(dict(
                detect_margin_median=float(np.nanmedian(frame["detect_margin"])),
                ident_margin_median=(float(np.nanmedian(frame["ident_margin"]))
                                     if estimable else float("nan")),
                jointly_evaluable_fraction=(float(frame["jointly_evaluable"].mean())
                                            if estimable else float("nan")),
                frac_R_below_one=fraction_below_one(frame["ratio_R"]),
                ratio_R_median=float(np.nanmedian(frame["ratio_R"])),
                n_excluded=len(report.excluded),
                # What each gate actually used. resolution_profiles clamps the detection
                # window to min(n_seeds, 8) while the competition side uses every split, so a
                # single recorded n_seeds would describe only one of the two axes.
                n_seeds_competition=_args.n_seeds,
                n_seeds_detection=int(window.n_seeds),
                competition_estimable=estimable,
                seed=_args.seed,
                control_label=control,
                seconds=round(time.time() - t0, 1),
            ))
            summary_rows.append(summary)
            print(f"[{gene:6} {space:3} depth {depth:3}] n={n:4} "
                  f"detection={report.detection_fraction:.3f} "
                  f"identification={report.identification_fraction:.3f} "
                  f"jointly={summary['jointly_evaluable_fraction']:.3f} "
                  f"R<1={summary['frac_R_below_one']:.3f} "
                  f"[{summary['seconds']:.1f}s]", flush=True)

            checkpoint(complete=False)

        del X

checkpoint(complete=True)

print(f"\n{sum(len(f) for f in unit_rows)} unit rows, {len(summary_rows)} summary rows "
      f"-> {OUT_DIR}")
print(f"total {round(time.time() - started, 1)} s")
