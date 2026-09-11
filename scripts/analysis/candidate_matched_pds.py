#!/usr/bin/env python3
"""Arm B of ``docs/PREREG_CANDIDATE_AND_POWER_v1.md``: does the JAK1 gap survive when the
candidate set is matched in size and in local geometry?

The four datasets rank against pools of very different size (JAK1 26, KRAS 93, TP53 98,
GATA1 255), and JAK1, the dataset carrying the paper's model-limited claim, has the
smallest. The objection is that its reproducibility reference is high because it has few
competitors. This resamples the candidate pool and scores every arm through the harness
the paper already uses, so the reference and the models move together and the reported
quantity is the gap.

Two sub-arms, answering different objections:

``random``
    A uniform subsample of the same-gene pool, always retaining the true target. Tests rank
    granularity and variance. Scored against the canonical constructions, so ``K = full``
    reproduces the committed per-variant scores exactly and that is a gate, not a hope.
    The subsample is drawn exactly rather than by enumerating candidates: PDS is a function
    of how many competitors fall closer than, tied with, and farther than the target, so
    the counts in a uniform draw are multivariate hypergeometric and are sampled directly.

``hard``
    The ``K - 1`` competitors nearest the target, which is the geometry Fig. 5g identified
    as governing. Selection uses the prediction half and evaluation the disjoint truth
    half, so for a **model** the competitors are chosen on cells the score never sees. The
    measurement arm is not clean here, because its prediction is the half the selection
    used; that shared draw biases its ``hard`` numbers by an amount and in a direction this
    analysis does not establish, so the arm is read for its shape and for the models rather
    than for the measurement's level. Every ``hard`` result is read against its own
    ``K = full``, not against the published value.

The decisive comparison is registered in advance and runs the direction that costs JAK1
nothing to fail: TP53, KRAS and GATA1 are down-matched to ``K = 25``, JAK1's pool size. If
they stay at chance there while JAK1 does not, pool size is not the explanation.

Usage:
    python scripts/analysis/candidate_matched_pds.py [--base BASE] --out OUT
    python scripts/analysis/candidate_matched_pds.py --out OUT --genes JAK1 --n-seed 3

    --base     directory holding the per-gene arrays (env: PERTRESOLVE_DATA).
    --out      directory receiving the two tables; may not be inside results/.
    --genes    restrict to these genes, comma separated (default: all four).
    --preds    directory of per-method prediction .npz (default: <base>/preds5).
    --n-seed   cell-subsample seeds (default 15, the canonical value).
    --n-draw   candidate draws per seed in the random arm (default 50).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import (  # noqa: E402
    add_harness_to_path,
    reject_repo_results,
    require_inputs,
    resolve_base,
)

NSUB = 300
NBOOT = 2000
TOL = 1e-12
WT_TAGS = ("WT", "wt", "WT_control")

#: Candidate-set sizes. 25 is JAK1's pool and the size the registered down-match uses;
#: sizes larger than a gene's pool are skipped rather than clipped, so a row always means
#: what its K says.
K_LADDER = (5, 10, 20, 25, 45, 90)

MEASUREMENT = "split-half measurement"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=None,
                    help="directory holding the per-gene arrays (env: PERTRESOLVE_DATA)")
    ap.add_argument("--out", required=True, metavar="OUT_DIR",
                    help="directory receiving candidate_matched_pds.csv.gz and "
                         "candidate_matched_summary.csv; may not be inside results/")
    ap.add_argument("--genes", default=None,
                    help="comma-separated subset of genes (default: all four)")
    ap.add_argument("--preds", default=None,
                    help="directory of per-method prediction .npz (default: <base>/preds5)")
    ap.add_argument("--n-seed", type=int, default=15,
                    help="cell-subsample seeds (default: %(default)s)")
    ap.add_argument("--n-draw", type=int, default=50,
                    help="candidate draws per seed, random arm (default: %(default)s)")
    return ap.parse_args()


def norm(M: np.ndarray) -> np.ndarray:
    """Row-normalise, leaving a zero row at zero rather than dividing by zero."""
    return M / (np.linalg.norm(M, axis=1, keepdims=True) + 1e-12)


def real_deltas_seed(H, gene_cells, genes, seed: int) -> dict:
    """Canonical per-variant pseudobulk profiles, verbatim from score_definitive.py."""
    rd = {}
    for g in genes:
        X, lab = gene_cells[g]
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


def halves_delta(H, gene_cells, genes, seed: int) -> tuple[dict, dict]:
    """Disjoint split-half deltas, verbatim from oracle_ceiling.py.

    truth[v] from one half plus a wild-type half, pred[v] from the disjoint half plus the
    complementary wild-type half. The wild-type split is drawn once per (gene, seed) and
    shared across variants; drawing it per variant collapses PDS far below chance.
    """
    truth, pred = {}, {}
    for g in genes:
        X, lab = gene_cells[g]
        rng = np.random.RandomState(1000 + 7 * seed + H._GENE_SEED[g])
        wt = np.where(np.isin(lab, WT_TAGS))[0]
        if len(wt) == 0:
            wm_t = wm_p = np.zeros(X.shape[1], np.float32)
        else:
            rng.shuffle(wt)
            h = len(wt) // 2
            wm_t = X[wt[:h][:NSUB]].mean(0)
            wm_p = X[wt[h:][:NSUB]].mean(0)
        for v in np.unique(lab):
            if v in WT_TAGS:
                continue
            idx = np.where(lab == v)[0]
            if len(idx) < 5:
                continue
            rng.shuffle(idx)
            h = len(idx) // 2
            truth[f"{g}__{v}"] = (X[idx[:h][:NSUB]].mean(0) - wm_t).astype(np.float32)
            pred[f"{g}__{v}"] = (X[idx[h:][:NSUB]].mean(0) - wm_p).astype(np.float32)
    return truth, pred


def pds_full(drow: np.ndarray, ci: int) -> float:
    """Tie-aware mid-rank PDS against the whole row, verbatim from score_definitive.py."""
    n = len(drow)
    td = drow[ci]
    if not np.isfinite(td):
        return 0.5
    less = int((drow < td - TOL).sum())
    eq = int((np.abs(drow - td) <= TOL).sum())
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1) if n > 1 else 0.5


def pds_random_k(drow: np.ndarray, ci: int, ks, rng, n_draw: int) -> dict[int, float]:
    """Mean PDS over uniform candidate subsamples of each size, target always retained.

    PDS depends on the competitors only through how many are closer, tied, and farther, so
    a uniform draw of ``K - 1`` of them gives multivariate hypergeometric counts. Sampling
    those counts directly is exact and avoids materialising a candidate list per draw.
    """
    td = drow[ci]
    if not np.isfinite(td):
        return {k: 0.5 for k in ks}
    others = np.delete(drow, ci)
    n_less = int((others < td - TOL).sum())
    n_tie = int((np.abs(others - td) <= TOL).sum())
    colors = np.array([n_less, n_tie, len(others) - n_less - n_tie])
    out = {}
    for k in ks:
        if k - 1 > colors.sum():
            continue
        drawn = rng.multivariate_hypergeometric(colors, k - 1, size=n_draw)
        out[k] = float(np.mean(1.0 - (drawn[:, 0] + drawn[:, 1] / 2.0) / (k - 1)))
    return out


def pds_hard_k(drow: np.ndarray, ci: int, sel: np.ndarray, ks) -> dict[int, float]:
    """PDS against the ``K - 1`` competitors nearest the target under ``sel``.

    ``sel`` is a selection distance from the target to every candidate, measured on cells
    disjoint from the evaluation truth, so the choice of competitor and the score that
    competitor produces do not share sampling noise.
    """
    order = np.argsort(np.delete(sel, ci), kind="stable")
    others_idx = np.delete(np.arange(len(drow)), ci)
    out = {}
    for k in ks:
        if k - 1 > len(order):
            continue
        keep = np.concatenate(([ci], others_idx[order[:k - 1]]))
        out[k] = pds_full(drow[keep], 0)
    return out


def bootstrap(values: np.ndarray, rng) -> tuple[float, float]:
    if len(values) < 2:
        return (np.nan, np.nan)
    draws = rng.integers(0, len(values), size=(NBOOT, len(values)))
    means = values[draws].mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


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
    genes = args.genes.split(",") if args.genes else list(H.GENES)
    unknown = [g for g in genes if g not in H.GENES]
    if unknown:
        raise SystemExit(f"unknown genes {unknown}; known: {list(H.GENES)}")

    methods = {p.stem: np.load(p, allow_pickle=True)
               for p in sorted(preds_dir.glob("*.npz"))}
    gene_cells = {g: (lambda X, lab: (X, np.asarray(lab)))(*H.load_gene(g)) for g in genes}
    splits = {g: sorted({k.split("__")[1] for k in methods[next(iter(methods))].files
                         if k.startswith(f"{g}__")}) for g in genes}
    cand = {(g, s): H.split_vars(g, s)[0] + H.split_vars(g, s)[1]
            for g in genes for s in splits[g]}
    print(f"{len(methods)} prediction files, {len(genes)} genes, "
          f"{args.n_seed} seeds, {args.n_draw} candidate draws per seed")

    # (gene, split, arm, method, K, variant) -> list of per-seed means
    acc: dict = {}

    def record(gene, split, arm, method, k, variant, value):
        acc.setdefault((gene, split, arm, method, k, variant), []).append(value)

    for seed in range(args.n_seed):
        rng = np.random.default_rng(20260828 + seed)
        real = real_deltas_seed(H, gene_cells, genes, seed)
        truth_h, pred_h = halves_delta(H, gene_cells, genes, seed)
        for g in genes:
            for s in splits[g]:
                cv = [v for v in cand[(g, s)] if f"{g}__{v}" in real]
                if len(cv) < 2:
                    continue
                idx = {v: i for i, v in enumerate(cv)}
                ks = [k for k in K_LADDER if k <= len(cv)]
                te_all = [v for v in H.split_vars(g, s)[1] if v in idx]

                # ---- measurement arm, canonical halves construction ----
                Tm_h = norm(np.stack([truth_h[f"{g}__{u}"] for u in cv]))
                Pm_h = norm(np.stack([pred_h[f"{g}__{v}"] for v in te_all]))
                D_h = 1 - Pm_h @ Tm_h.T
                Sel = 1 - Pm_h @ norm(np.stack([pred_h[f"{g}__{u}"]
                                                for u in cv])).T
                for i, v in enumerate(te_all):
                    ci = idx[v]
                    record(g, s, "random", MEASUREMENT, len(cv), v, pds_full(D_h[i], ci))
                    record(g, s, "hard", MEASUREMENT, len(cv), v, pds_full(D_h[i], ci))
                    for k, val in pds_random_k(D_h[i], ci, ks, rng, args.n_draw).items():
                        record(g, s, "random", MEASUREMENT, k, v, val)
                    for k, val in pds_hard_k(D_h[i], ci, Sel[i], ks).items():
                        record(g, s, "hard", MEASUREMENT, k, v, val)

                # ---- model arm ----
                Tm_r = norm(np.stack([real[f"{g}__{u}"] for u in cv]))
                for name, store in methods.items():
                    te = [v for v in te_all if f"{g}__{s}__{v}" in store.files]
                    if not te:
                        continue
                    raw = np.stack([np.asarray(store[f"{g}__{s}__{v}"], dtype=float)
                                    for v in te])
                    finite = np.all(np.isfinite(raw), axis=1)
                    if not finite.any():
                        continue          # refused; coverage is reported by arm C
                    te = [v for v, ok in zip(te, finite) if ok]
                    Pm = norm(raw[finite])
                    D_r = 1 - Pm @ Tm_r.T          # canonical truth: reproduces the table
                    D_hm = 1 - Pm @ Tm_h.T         # half truth: the hard arm's baseline
                    for i, v in enumerate(te):
                        ci = idx[v]
                        record(g, s, "random", name, len(cv), v, pds_full(D_r[i], ci))
                        record(g, s, "hard", name, len(cv), v, pds_full(D_hm[i], ci))
                        for k, val in pds_random_k(D_r[i], ci, ks, rng,
                                                   args.n_draw).items():
                            record(g, s, "random", name, k, v, val)
                        # Selection depends on the target variant only, never on the
                        # method, so the row is reused from the measurement arm rather
                        # than rebuilding the candidate matrix once per method.
                        sel_row = Sel[te_all.index(v)]
                        for k, val in pds_hard_k(D_hm[i], ci, sel_row, ks).items():
                            record(g, s, "hard", name, k, v, val)
        print(f"seed {seed} done", flush=True)

    rows = [dict(gene=g, split=s, arm=a, method=m, K=k, variant=v,
                 pds=float(np.mean(vals)), n_seed=len(vals))
            for (g, s, a, m, k, v), vals in acc.items()]
    per_var = pd.DataFrame(rows)
    per_var.to_csv(out_dir / "candidate_matched_pds.csv.gz", index=False,
                   float_format="%.6f", compression="gzip")

    summary = summarise(per_var, np.random.default_rng(0))
    summary.to_csv(out_dir / "candidate_matched_summary.csv", index=False)
    print(f"\nwrote {out_dir / 'candidate_matched_pds.csv.gz'} ({len(per_var):,} rows) "
          f"and {out_dir / 'candidate_matched_summary.csv'} ({len(summary):,} rows)")

    check_gates(per_var, genes)


def summarise(per_var: pd.DataFrame, rng) -> pd.DataFrame:
    """Aggregate over distinct variants, the unit the manuscript bootstraps.

    A variant can be held out in several splits, so it is averaged across them before
    entering the bootstrap rather than contributing once per split.
    """
    out = []
    keys = ["gene", "arm", "method", "K"]
    per_v = (per_var.groupby(keys + ["variant"], as_index=False).pds.mean())
    for key, grp in per_v.groupby(keys, sort=True):
        vals = grp.pds.to_numpy()
        lo, hi = bootstrap(vals, rng)
        out.append(dict(zip(keys, key)) | dict(
            n_var=len(vals), pds=round(float(vals.mean()), 5),
            ci_lo=round(lo, 5) if np.isfinite(lo) else np.nan,
            ci_hi=round(hi, 5) if np.isfinite(hi) else np.nan))
    return pd.DataFrame(out)


def check_gates(per_var: pd.DataFrame, genes: list[str]) -> None:
    """At the full pool the random arm must reproduce the published reference.

    The gate is the split-half reference in ``results/canonical/oracle_ceiling.csv``,
    aggregated the way ``oracle_ceiling.summarize`` does: seed-averaged per variant, meaned
    within each (gene, split) cell, then meaned over cells. Getting that aggregation wrong
    is itself a way to fail silently, so ``n_var_scored`` is checked as well as the value.

    The per-variant model table is compared but does **not** gate. It is not row
    reproducible from ``score_definitive.py``: re-running that generator, both here and on
    the machine that produced it, returns values that agree with the committed table in
    aggregate and not row by row. See docs/RESULT_CANDIDATE_AND_POWER_v1.md. The comparison
    is printed so the size of that disagreement travels with every run.
    """
    print("\n=== gate: K = full split-half reference vs oracle_ceiling.csv ===")
    repo = Path(__file__).resolve().parents[2]
    ok = True

    full = per_var[per_var.arm == "random"].copy()
    full = full[full.K == full.groupby(["gene", "split"]).K.transform("max")]

    ceiling = repo / "results" / "canonical" / "oracle_ceiling.csv"
    if ceiling.exists():
        ref = pd.read_csv(ceiling).set_index("scope")
        meas = full[full.method == MEASUREMENT]
        cell = meas.groupby(["gene", "split"], as_index=False).pds.mean()
        got = cell.groupby("gene").pds.mean()
        n_scored = meas.groupby("gene").size()
        for g in genes:
            if g not in ref.index or g not in got.index:
                continue
            d = abs(float(got[g]) - float(ref.loc[g, "PDS_oracle"]))
            n_ok = int(n_scored[g]) == int(ref.loc[g, "n_var_scored"])
            passed = d <= 0.005 and n_ok
            ok &= passed
            print(f"{'PASS' if passed else 'FAIL'}  {g} {got[g]:.4f} vs published "
                  f"{float(ref.loc[g, 'PDS_oracle']):.4f}; scored {int(n_scored[g])} vs "
                  f"{int(ref.loc[g, 'n_var_scored'])}")

    print("\n=== diagnostic, not a gate: model scores vs per_seed_variant_pds.csv.gz ===")
    committed = repo / "results" / "canonical" / "per_seed_variant_pds.csv.gz"
    if committed.exists():
        pub = (pd.read_csv(committed)
                 .groupby(["method", "gene", "split", "variant"], as_index=False).pds.mean()
                 .rename(columns={"pds": "pds_pub"}))
        merged = full[full.method != MEASUREMENT].merge(
            pub, on=["method", "gene", "split", "variant"], how="inner")
        if merged.empty:
            print("no overlap with the committed per-variant scores")
        else:
            d = (merged.pds - merged.pds_pub).abs()
            mine = (merged.groupby(["method", "gene", "split"], as_index=False).pds.mean()
                          .groupby("method").pds.mean().mean())
            theirs = (merged.groupby(["method", "gene", "split"], as_index=False)
                            .pds_pub.mean().groupby("method").pds_pub.mean().mean())
            print(f"{len(merged):,} matched rows: median |diff| {d.median():.5f}, "
                  f"{(d > 0.05).sum():,} rows over 0.05, max {d.max():.5f}")
            print(f"aggregate over all methods: recompute {mine:.4f} vs committed "
                  f"{theirs:.4f} (difference {abs(mine - theirs):.4f})")
    if not ok:
        raise SystemExit("the reference gate failed; the tables above are not to be read")


if __name__ == "__main__":
    main()
