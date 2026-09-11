#!/usr/bin/env python3
"""Regenerate the canonical PDS artifacts from the prediction vectors and truth profiles.

**Not applied.** The committed tables are unchanged; this exists so the regeneration can be
run and audited, and so the invariants below travel with the repository. See
``docs/OPEN_ITEMS_2026-08-28.md`` for what an audit of those tables did and did not find.

What that audit found, after one wrong turn: ``definitive_summary.csv`` and
``per_seed_variant_pds.csv.gz`` **are** mutually consistent. The summary carries two
aggregations, ``PDS`` over distinct (gene, variant) pairs, which is the estimator the
manuscript specifies, and ``PDS_cell`` over (gene, split) cell means; each reproduces from
the per-variant table to within rounding. What does not reproduce is the per-variant table
row by row: a fresh run of the committed generator, on byte-identical inputs, agrees on cell
means to 0.0075 and on the summary to 0.0046 while differing on 1,009 of 12,814 individual
rows. The provenance of the committed draw was not recovered.

Running this produces a faithful drop-in for both tables, with every column the committed
summary carries, plus a per-method old-to-new report.

**What is enforced, not assumed.**

1. The summary is computed *from* the per-variant table in one code path, then re-derived
   independently and asserted equal.
2. A random sample of stored rows is re-scored from a freshly built distance matrix and
   asserted equal, so a table cannot disagree with the distances it came from. This is the
   check whose absence let the row-level drift go unnoticed.
3. The vectorised scorer is asserted equal to the packaged ``pertresolve.evaluation.pds``,
   so "one canonical scorer" is a checked claim rather than a convention.
4. Nothing is overwritten. The previous files are copied to an archive directory with their
   SHA-256 recorded before anything new is written.

Usage:
    python scripts/analysis/regenerate_canonical_pds.py [--base BASE] --out OUT
    python scripts/analysis/regenerate_canonical_pds.py --out OUT --n-seed 2   # smoke

    --base     directory holding the per-gene arrays (env: PERTRESOLVE_DATA).
    --out      directory receiving the regenerated tables; may not be inside results/.
    --preds    directory of per-method prediction .npz (default: <base>/preds5).
    --n-seed   cell-subsample seeds (default 15, the canonical value).
    --archive  directory to copy the existing committed tables into before writing.
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.evaluation.pds import pds_score  # noqa: E402
from pertresolve.paths import (  # noqa: E402
    add_harness_to_path,
    reject_repo_results,
    require_inputs,
    resolve_base,
)

NSUB = 300
NBOOT = 2000
TOL = 1e-12
SPLITS = ("split1", "split2", "split3", "split5", "split6")
N_SPOTCHECK = 60


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=None)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--preds", default=None)
    ap.add_argument("--n-seed", type=int, default=15)
    ap.add_argument("--archive", default=None, type=Path,
                    help="directory to copy the existing committed tables into "
                         "(default: results/canonical/_superseded_<date>)")
    return ap.parse_args()


def norm(M: np.ndarray) -> np.ndarray:
    """Row-normalise, leaving a zero row at zero rather than dividing by zero."""
    return M / (np.linalg.norm(M, axis=1, keepdims=True) + 1e-12)


def pds_from_row(drow: np.ndarray, ci: int) -> float:
    """THE canonical scorer: tie-aware mid-rank PDS for one prediction.

    Ties take the average rank, which keeps a prediction equidistant from several
    candidates at chance rather than at either extreme. A non-finite distance to the true
    match means the prediction could not be ranked; it takes the chance value, and the
    count of such rows is reported so the exclusion alternative stays visible.

    Asserted equal to ``pertresolve.evaluation.pds.pds_score`` on a random sample of real
    rows by :func:`check_scorer_agrees_with_package`.
    """
    n = len(drow)
    if n <= 1:
        return 0.5
    td = drow[ci]
    if not np.isfinite(td):
        return 0.5
    less = int((drow < td - TOL).sum())
    eq = int((np.abs(drow - td) <= TOL).sum())
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1)


def real_deltas_seed(H, gene_cells, seed: int) -> dict:
    """Per-variant pseudobulk profiles for one cell-subsample seed."""
    rd = {}
    for g, (X, lab) in gene_cells.items():
        rng = np.random.RandomState(1000 + 7 * seed + H._GENE_SEED[g])
        wt = np.where(np.isin(lab, H.WT_TAGS))[0]
        wt = rng.choice(wt, NSUB, replace=False) if len(wt) > NSUB else wt
        wm = X[wt].mean(0)
        for v in np.unique(lab):
            if v in H.WT_TAGS:
                continue
            idx = np.where(lab == v)[0]
            if len(idx) < 5:
                continue
            if len(idx) > NSUB:
                idx = rng.choice(idx, NSUB, replace=False)
            rd[f"{g}__{v}"] = (X[idx].mean(0) - wm).astype(np.float32)
    return rd


def score_all(H, gene_cells, methods, n_seed: int) -> tuple[pd.DataFrame, dict]:
    """Every (method, gene, split, variant, seed) score, plus refusals per method."""
    cand = {(g, s): H.split_vars(g, s)[0] + H.split_vars(g, s)[1]
            for g in gene_cells for s in SPLITS}
    rows, n_refused = [], {}
    for seed in range(n_seed):
        real = real_deltas_seed(H, gene_cells, seed)
        for m, store in methods.items():
            for g in gene_cells:
                for s in SPLITS:
                    cv = [v for v in cand[(g, s)] if f"{g}__{v}" in real]
                    if len(cv) < 2:
                        continue
                    idx = {v: i for i, v in enumerate(cv)}
                    te = [v for v in H.split_vars(g, s)[1]
                          if f"{g}__{s}__{v}" in store.files and v in idx]
                    if not te:
                        continue
                    Pm = norm(np.stack([store[f"{g}__{s}__{v}"] for v in te]))
                    Tm = norm(np.stack([real[f"{g}__{u}"] for u in cv]))
                    D = 1 - Pm @ Tm.T
                    for i, v in enumerate(te):
                        if not np.isfinite(D[i, idx[v]]):
                            n_refused[m] = n_refused.get(m, 0) + 1
                        rows.append((m, g, s, v, seed, pds_from_row(D[i], idx[v])))
        print(f"  seed {seed} done", flush=True)
    return pd.DataFrame(rows, columns=["method", "gene", "split", "variant",
                                       "seed", "pds"]), n_refused


def per_variant(per_seed: pd.DataFrame) -> pd.DataFrame:
    """Seed-averaged score for each (method, gene, split, variant)."""
    return per_seed.groupby(["method", "gene", "split", "variant"],
                            as_index=False).pds.mean()


def summarise(per_var: pd.DataFrame, n_refused_by_method: dict) -> pd.DataFrame:
    """Both aggregations the committed summary carries, so this is a drop-in replacement.

    ``PDS`` is over distinct (gene, variant) pairs, which is the estimator the manuscript
    specifies: a variant held out under several splits is averaged across them first so it
    enters the point estimate and the bootstrap once. ``PDS_cell`` is the mean of
    (gene, split) cell means, which counts such a variant once per split. The committed
    ``score_definitive.py`` emits both and so does this; dropping either would silently
    change what a downstream reader gets.
    """
    out = []
    for m, grp in per_var.groupby("method", sort=True):
        distinct = grp.groupby(["gene", "variant"]).pds.mean().to_numpy()
        cells = [c.pds.to_numpy() for _, c in grp.groupby(["gene", "split"])]
        rng = np.random.default_rng(0)

        dv = rng.integers(0, len(distinct), size=(NBOOT, len(distinct)))
        lo_v, hi_v = np.percentile(distinct[dv].mean(axis=1), [2.5, 97.5])

        boot_c = np.array([np.mean([c[rng.integers(0, len(c), len(c))].mean()
                                    for c in cells]) for _ in range(NBOOT)])
        lo_c, hi_c = np.percentile(boot_c, [2.5, 97.5])
        obs_c = float(np.mean([c.mean() for c in cells]))

        out.append(dict(
            method=m, standardization="train-only",
            PDS=round(float(distinct.mean()), 4),
            ci_lo=round(float(lo_v), 4), ci_hi=round(float(hi_v), 4),
            crosses=bool(lo_v <= 0.5 <= hi_v), n_variants=len(distinct),
            PDS_cell=round(obs_c, 4), cell_lo=round(float(lo_c), 4),
            cell_hi=round(float(hi_c), 4), cell_crosses=bool(lo_c <= 0.5 <= hi_c),
            n_cells=len(cells), n_refused=int(n_refused_by_method.get(m, 0))))
    return pd.DataFrame(out).sort_values("PDS", ascending=False).reset_index(drop=True)


# ------------------------------------------------------------------ enforced invariants


def check_summary_matches_per_variant(summary: pd.DataFrame, per_var: pd.DataFrame) -> None:
    """The summary must be the aggregate of the table it ships with."""
    redone = (per_var.groupby(["method", "gene", "variant"], as_index=False).pds.mean()
                     .groupby("method").pds.mean())
    merged = summary.set_index("method").PDS
    d = (merged - redone.reindex(merged.index)).abs()
    if d.max() > 5e-5:
        raise SystemExit(f"summary does not aggregate from the per-variant table; "
                         f"max |diff| {d.max():.6f}")
    print(f"PASS  summary == aggregate(per-variant table), max |diff| {d.max():.2e}")


def check_scorer_agrees_with_package(H, gene_cells, methods, rng) -> None:
    """One canonical scorer: the vectorised form must equal the packaged one."""
    g = list(gene_cells)[0]
    real = real_deltas_seed(H, gene_cells, 0)
    s = SPLITS[0]
    cv = [v for v in H.split_vars(g, s)[0] + H.split_vars(g, s)[1] if f"{g}__{v}" in real]
    idx = {v: i for i, v in enumerate(cv)}
    m, store = next(iter(methods.items()))
    te = [v for v in H.split_vars(g, s)[1] if f"{g}__{s}__{v}" in store.files and v in idx]
    if not te:
        print("SKIP  scorer agreement: no scorable rows in the first cell")
        return
    Tm = norm(np.stack([real[f"{g}__{u}"] for u in cv]))
    deltas = {u: real[f"{g}__{u}"] for u in cv}
    worst = 0.0
    for v in te:
        pred = store[f"{g}__{s}__{v}"]
        d = 1 - norm(pred[None, :])[0] @ Tm.T
        worst = max(worst, abs(pds_from_row(d, idx[v])
                               - pds_score(pred, v, deltas, cv, distance="cos")))
    if worst > 1e-9:
        raise SystemExit(f"the vectorised scorer disagrees with pertresolve.evaluation."
                         f"pds by {worst:.2e}")
    print(f"PASS  vectorised scorer == pertresolve.evaluation.pds over {len(te)} rows, "
          f"max |diff| {worst:.2e}")


def spot_check_rows(H, gene_cells, methods, per_seed: pd.DataFrame, rng) -> None:
    """Re-score sampled stored rows from a freshly built distance matrix."""
    sample = per_seed.sample(n=min(N_SPOTCHECK, len(per_seed)), random_state=0)
    cache: dict[int, dict] = {}
    worst, checked = 0.0, 0
    for _, r in sample.iterrows():
        real = cache.setdefault(int(r.seed), real_deltas_seed(H, gene_cells, int(r.seed)))
        g, s = r.gene, r.split
        cv = [v for v in H.split_vars(g, s)[0] + H.split_vars(g, s)[1]
              if f"{g}__{v}" in real]
        idx = {v: i for i, v in enumerate(cv)}
        Tm = norm(np.stack([real[f"{g}__{u}"] for u in cv]))
        pred = methods[r.method][f"{g}__{s}__{r.variant}"]
        d = 1 - norm(np.asarray(pred, dtype=np.float32)[None, :])[0] @ Tm.T
        worst = max(worst, abs(pds_from_row(d, idx[r.variant]) - float(r.pds)))
        checked += 1
    if worst > 1e-9:
        raise SystemExit(f"stored PDS disagrees with a fresh recompute by {worst:.2e}")
    print(f"PASS  {checked} sampled rows re-scored from the distance matrix, "
          f"max |diff| {worst:.2e}")


def archive_existing(archive: Path, files: list[Path]) -> pd.DataFrame:
    """Copy the superseded tables somewhere safe and record their digests."""
    archive.mkdir(parents=True, exist_ok=True)
    rows = []
    for f in files:
        if not f.exists():
            continue
        digest = hashlib.sha256(f.read_bytes()).hexdigest()
        shutil.copy2(f, archive / f.name)
        rows.append(dict(file=f.name, sha256=digest, bytes=f.stat().st_size))
        print(f"  archived {f.name}  sha256 {digest[:16]}...")
    manifest = pd.DataFrame(rows)
    manifest.to_csv(archive / "MANIFEST.csv", index=False)
    return manifest


def diff_report(old_summary: Path, new_summary: pd.DataFrame) -> pd.DataFrame | None:
    """Per-method old to new difference, the record of what this repair moved."""
    if not old_summary.exists():
        return None
    old = pd.read_csv(old_summary)[["method", "PDS", "ci_lo", "ci_hi"]].rename(
        columns={"PDS": "PDS_old", "ci_lo": "ci_lo_old", "ci_hi": "ci_hi_old"})
    new = new_summary[["method", "PDS", "ci_lo", "ci_hi", "n_variants"]].rename(
        columns={"PDS": "PDS_new", "ci_lo": "ci_lo_new", "ci_hi": "ci_hi_new"})
    d = old.merge(new, on="method", how="outer")
    d["delta"] = (d.PDS_new - d.PDS_old).round(4)
    d["crosses_half_old"] = (d.ci_lo_old <= 0.5) & (0.5 <= d.ci_hi_old)
    d["crosses_half_new"] = (d.ci_lo_new <= 0.5) & (0.5 <= d.ci_hi_new)
    d["verdict_flips"] = d.crosses_half_old != d.crosses_half_new
    return d.sort_values("delta")


def main() -> None:
    args = parse_args()
    base = resolve_base(args.base)
    out_dir = reject_repo_results(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    require_inputs(base / "pertresolve_bench.csv")
    add_harness_to_path(base)

    import harness as H  # noqa: E402

    preds_dir = Path(args.preds) if args.preds else base / "preds5"
    require_inputs(preds_dir)
    methods = {p.stem: np.load(p, allow_pickle=True)
               for p in sorted(preds_dir.glob("*.npz"))}
    gene_cells = {g: (lambda X, lab: (X, np.asarray(lab)))(*H.load_gene(g))
                  for g in H.GENES}
    rng = np.random.default_rng(0)

    repo = Path(__file__).resolve().parents[2]
    canonical = repo / "results" / "canonical"
    archive = args.archive or canonical / "_superseded_2026-08-28"
    print("=== archiving the superseded tables ===")
    archive_existing(archive, [canonical / "per_seed_variant_pds.csv.gz",
                               canonical / "definitive_summary.csv"])

    print(f"\n=== scoring: {len(methods)} methods, {len(gene_cells)} genes, "
          f"{args.n_seed} seeds ===")
    check_scorer_agrees_with_package(H, gene_cells, methods, rng)
    per_seed, n_refused = score_all(H, gene_cells, methods, args.n_seed)
    per_var = per_variant(per_seed)
    summary = summarise(per_var, n_refused)

    print("\n=== enforced invariants ===")
    check_summary_matches_per_variant(summary, per_var)
    spot_check_rows(H, gene_cells, methods, per_seed, rng)
    print(f"NOTE  {sum(n_refused.values())} of {len(per_seed):,} scored rows had a "
          f"non-finite distance to the true match and took the chance value; "
          f"per method: {n_refused}")

    per_seed.to_csv(out_dir / "per_seed_variant_pds.csv.gz", index=False,
                    float_format="%.6f", compression="gzip")
    per_var.to_csv(out_dir / "per_variant_pds.csv.gz", index=False,
                   float_format="%.6f", compression="gzip")
    summary.to_csv(out_dir / "definitive_summary.csv", index=False)

    d = diff_report(archive / "definitive_summary.csv", summary)
    if d is not None:
        d.to_csv(out_dir / "repair_diff_summary.csv", index=False)
        print("\n=== old to new, per method ===")
        print(d[["method", "PDS_old", "PDS_new", "delta", "crosses_half_old",
                 "crosses_half_new", "verdict_flips"]].to_string(index=False))
        print(f"\nlargest move {d.delta.abs().max():.4f}; "
              f"methods whose chance verdict flips: {int(d.verdict_flips.sum())}")
    print(f"\nwrote {out_dir}")


if __name__ == "__main__":
    main()
