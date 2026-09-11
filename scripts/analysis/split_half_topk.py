#!/usr/bin/env python3
"""Split-half Top-k recovery: does a disjoint measurement put the true allele in a shortlist?

This is the same experiment as ``oracle_ceiling.py``, read differently. That script turns
each held-out variant's rank into a percentile (PDS) and reports the mean; this one asks
the question an experimenter actually faces, which is whether the correct allele reaches a
shortlist of the k most similar candidates.

**No new similarity, no new candidate set, no new inclusion rule.** The build-half response,
the evaluation-half candidate profiles, the L2-normalised cosine distance, the tie handling,
the 15 seeds, the five splits, the 300-cell subsample cap and the five-cell minimum are all
taken verbatim from ``oracle_ceiling.py``. The only addition is that the rank the percentile
was computed from is retained instead of being discarded, and the same aggregation is applied
to an indicator instead of to a percentile.

Why it cannot be derived from the committed table. ``oracle_ceiling_per_variant.csv`` stores
the mean of 15 per-seed percentiles. Hit@k is the mean of 15 per-seed *indicators*, and the
mean of an indicator is not a function of the mean of its argument: a variant at rank 1 in
seven seeds and rank 60 in eight has a mean percentile near the middle of the pool and a
Hit@1 of 7/15, while a variant at rank 30 in all fifteen has a similar mean percentile and a
Hit@1 of 0. Inverting the committed mean yields a mean rank, which is a different quantity
and would answer a different question. Hence this re-run.

**Ties.** ``pds_row`` averages the ranks of tied candidates. The matching convention for an
indicator is the expected hit under uniform tie-breaking: with ``less`` candidates strictly
closer and ``eq`` candidates at the target's distance (the target included), the rank is
uniform on ``{less+1, ..., less+eq}``, so ``P(rank <= k) = clip(k - less, 0, eq) / eq``.
Without ties this reduces to ``1(less < k)``, the plain indicator. On the committed run one
variant produced ties, GATA1 S274G, which appears in two splits and so in two unit rows.

**Aggregation** is ``summarize_variant``'s, the distinct-variant unit Methods declares as
primary and the one Fig. 3b plots: mean over seeds, then over splits within a variant, then
over distinct variants. The bootstrap resamples distinct variants, 2,000 draws from
``default_rng(0)``, the same construction as there.

The middle step is a formality on these data and is kept only so the two tables aggregate
identically. Each split's candidate set is ``train + test``, which for every split is the
gene's whole scored pool, and neither ``truth`` nor ``pred`` depends on the split, so a
variant held out under several splits receives the same score under each. The 518 rows of
the unit table therefore carry 308 distinct values, and the effective sample is the 308
distinct variants, not 518 independent units. The same is true of ``oracle_ceiling.py``'s
own per-variant table, which is why the two agree row for row.

**Null.** Under a uniform random ordering of a pool of ``N`` candidates the true allele
reaches the top k with probability ``min(k/N, 1)``. The pool differs by gene (98, 92, 254 and
26 candidates), so the null differs by gene and a single chance line would be wrong for three
of the four. It is computed per scored unit and averaged with the observed, then checked
against the closed form.

The script re-derives the committed ``pds`` column as it goes and refuses to write anything
if it disagrees, so the two tables are provably the same run of the same harness.

Run (compute box, not this repository's results tree):

    python scripts/analysis/split_half_topk.py \
        --base /path/to/processed-data --out /tmp/topk
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import (
    add_harness_to_path,
    reject_repo_results,
    repo_root,
    resolve_base,
)

#: Verbatim from oracle_ceiling.py. Changing any of these makes this a different harness.
NSUB = 300
NSEED = 15
NBOOT = 2000
SPLITS = ["split1", "split2", "split3", "split5", "split6"]
WT_TAGS = ("WT", "wt", "WT_control")
MIN_CELLS = 5

#: The four shortlist lengths. Fixed here rather than exposed as a flag: a k chosen after
#: seeing the curves would be a threshold moved to fit the result.
KS = (1, 3, 5, 10)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=None,
                    help="directory holding the shared scorer and the gene arrays "
                         "(env: PERTRESOLVE_DATA)")
    ap.add_argument("--out", required=True, metavar="OUT_DIR",
                    help="directory receiving split_half_topk.csv; may not be inside the "
                         "repository's results/")
    return ap.parse_args()


def halves_delta(H, gene_cells, seed):
    """Disjoint split-half pseudobulk deltas. Copied from oracle_ceiling.py unchanged."""
    truth, pred = {}, {}
    for g in H.GENES:
        X, lab = gene_cells[g]
        rng = np.random.RandomState(1000 + 7 * seed + H._GENE_SEED[g])
        wt = np.where(np.isin(lab, WT_TAGS))[0]
        if len(wt) == 0:
            wm_t = wm_p = np.zeros(X.shape[1], np.float32)
        else:
            rng.shuffle(wt)
            h = len(wt) // 2
            wt_t, wt_p = wt[:h], wt[h:]
            if len(wt_t) > NSUB:
                wt_t = wt_t[:NSUB]
            if len(wt_p) > NSUB:
                wt_p = wt_p[:NSUB]
            wm_t = X[wt_t].mean(0)
            wm_p = X[wt_p].mean(0)
        for v in np.unique(lab):
            if v in WT_TAGS:
                continue
            idx = np.where(lab == v)[0]
            if len(idx) < MIN_CELLS:
                continue
            rng.shuffle(idx)
            h = len(idx) // 2
            it, ip = idx[:h], idx[h:]
            if len(it) > NSUB:
                it = it[:NSUB]
            if len(ip) > NSUB:
                ip = ip[:NSUB]
            truth[f"{g}__{v}"] = (X[it].mean(0) - wm_t).astype(np.float32)
            pred[f"{g}__{v}"] = (X[ip].mean(0) - wm_p).astype(np.float32)
    return truth, pred


def norm(M):
    return M / (np.linalg.norm(M, axis=1, keepdims=True) + 1e-12)


def rank_counts(drow, ci):
    """Candidates strictly closer than the target, and candidates tied with it.

    ``eq`` counts the target itself, so it is at least 1. Both are integers, so the rank
    and the percentile below are exact rather than reconstructed from a float.
    """
    td = drow[ci]
    if not np.isfinite(td):
        return None
    less = int((drow < td - 1e-12).sum())
    eq = int((np.abs(drow - td) <= 1e-12).sum())
    return less, eq


def pds_from_counts(less, eq, n):
    """``oracle_ceiling.pds_row``, written in terms of the retained counts."""
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1) if n > 1 else 0.5


def hit_from_counts(less, eq, k):
    """Expected hit at k under uniform tie-breaking; the plain indicator when eq == 1."""
    return float(min(max(k - less, 0), eq)) / eq


def collect(H, gene_cells, cand):
    """One record per (gene, split, variant, seed): the counts the rank is made of."""
    acc: dict = {}
    pools: dict = {}
    non_finite = 0
    for seed in range(NSEED):
        truth, pred = halves_delta(H, gene_cells, seed)
        for g in H.GENES:
            for s in SPLITS:
                cv = [v for v in cand[(g, s)] if f"{g}__{v}" in truth]
                if len(cv) < 2:
                    continue
                idx = {v: i for i, v in enumerate(cv)}
                te = [v for v in H.split_vars(g, s)[1]
                      if f"{g}__{v}" in truth and v in idx]
                if not te:
                    continue
                pools.setdefault((g, s), set()).add(len(cv))
                Tm = norm(np.stack([truth[f"{g}__{u}"] for u in cv]))
                Pm = norm(np.stack([pred[f"{g}__{v}"] for v in te]))
                D = 1 - Pm @ Tm.T
                for i, v in enumerate(te):
                    counts = rank_counts(D[i], idx[v])
                    if counts is None:
                        non_finite += 1
                        continue
                    acc.setdefault((g, s), {}).setdefault(v, []).append(
                        (counts[0], counts[1], len(cv)))
    if non_finite:
        raise SystemExit(
            f"{non_finite} scored units produced a non-finite distance to the target. "
            "oracle_ceiling.py silently scores those as chance; a shortlist indicator has "
            "no such convention, so this run stops rather than inventing one.")
    return acc, pools


def per_unit_frame(acc) -> pd.DataFrame:
    """(gene, split, variant) rows carrying the seed-mean of every derived quantity."""
    rows = []
    for (g, s), by_variant in acc.items():
        for v, seeds in by_variant.items():
            n_cand = {c for _, _, c in seeds}
            if len(n_cand) != 1:
                raise SystemExit(
                    f"the candidate pool for {g}/{s}/{v} changed between seeds: {n_cand}. "
                    "The pool is defined by the split and the five-cell minimum, neither "
                    "of which depends on the seed, so this means the harness moved.")
            n = n_cand.pop()
            ranks = [less + (eq + 1) / 2.0 for less, eq, _ in seeds]
            row = dict(
                gene=g, split=s, variant=v, n_seeds=len(seeds),
                candidate_pool_size=n,
                mean_rank=float(np.mean(ranks)),
                median_rank=float(np.median(ranks)),
                pds=float(np.mean([pds_from_counts(less, eq, n) for less, eq, _ in seeds])),
            )
            for k in KS:
                row[f"hit_{k}"] = float(np.mean(
                    [hit_from_counts(less, eq, k) for less, eq, _ in seeds]))
                row[f"null_{k}"] = min(k / n, 1.0)
            rows.append(row)
    return pd.DataFrame(rows).sort_values(["gene", "split", "variant"]).reset_index(drop=True)


def summarize_variant(unit: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """oracle_ceiling.summarize_variant's aggregation, applied to each named column.

    Mean over splits within a distinct variant, then over distinct variants, with a
    2,000-draw bootstrap over distinct variants from ``default_rng(0)``.

    The interval this returns is a percentile bootstrap of the observed rate. It states how
    precisely the rate is known; it is not a test against the analytic chance rate, and it
    is only conservative for that purpose where it lies wholly above it, because its width
    scales with the observed rate. Fig. 3g quotes a fold value only in that direction.
    """
    out = []
    for scope in list(dict.fromkeys(unit["gene"])) + ["ALL"]:
        sub = unit if scope == "ALL" else unit[unit["gene"] == scope]
        distinct = sub.groupby(["gene", "variant"])[columns].mean()
        rng = np.random.default_rng(0)
        draw = rng.integers(0, len(distinct), size=(NBOOT, len(distinct)))
        record = dict(scope=scope, n_variants=len(distinct))
        for column in columns:
            values = distinct[column].to_numpy()
            lo, hi = np.percentile(values[draw].mean(axis=1), [2.5, 97.5])
            record[column] = float(values.mean())
            record[f"{column}_lo"] = float(lo)
            record[f"{column}_hi"] = float(hi)
        out.append(record)
    return pd.DataFrame(out)


def verify_against_committed(unit: pd.DataFrame) -> str:
    """Re-derived PDS must equal the committed split-half table, or nothing is written."""
    try:
        committed = repo_root() / "results" / "canonical" / "oracle_ceiling_per_variant.csv"
    except RuntimeError as exc:
        return f"skipped: the repository is not importable from here ({exc})"
    if not committed.exists():
        return f"skipped: {committed} is not present on this machine"
    ref = pd.read_csv(committed)
    key = ["gene", "split", "variant"]
    merged = ref.merge(unit[key + ["pds", "n_seeds"]], on=key, how="outer",
                       suffixes=("_ref", "_new"), indicator=True)
    unmatched = merged[merged["_merge"] != "both"]
    if len(unmatched):
        raise SystemExit(
            f"the re-run scored a different set of units than the committed table: "
            f"{len(unmatched)} rows do not match on (gene, split, variant).\n"
            f"{unmatched.head(12).to_string(index=False)}")
    delta = (merged["pds_ref"] - merged["pds_new"]).abs().max()
    if not (delta < 1e-9):
        worst = merged.assign(d=(merged["pds_ref"] - merged["pds_new"]).abs()) \
                      .nlargest(8, "d")[key + ["pds_ref", "pds_new", "d"]]
        raise SystemExit(
            "the re-derived PDS does not reproduce the committed split-half table "
            f"(max |difference| = {delta:.3e}), so this is not the same harness and the "
            f"Top-k numbers would not be comparable with Fig. 3b:\n{worst.to_string(index=False)}")
    return f"reproduced {len(merged)} committed per-variant PDS values, max |diff| = {delta:.2e}"


def check(unit: pd.DataFrame, summary: pd.DataFrame, table: pd.DataFrame) -> None:
    """Definitional and completeness checks only.

    Nothing here asserts which gene does well. The panel's reading has to come out of the
    numbers; these check that the numbers mean what the axis says they mean.
    """
    if sorted(table["k"].unique()) != sorted(KS):
        raise SystemExit(f"expected k in {KS}, found {sorted(table['k'].unique())}")
    bad = table[(table["observed_hit_rate"] < 0) | (table["observed_hit_rate"] > 1)]
    if len(bad):
        raise SystemExit(f"a hit rate left [0, 1]:\n{bad.to_string(index=False)}")
    for gene, rows in table.groupby("scope"):
        rows = rows.sort_values("k")
        observed = rows["observed_hit_rate"].to_numpy()
        if np.any(np.diff(observed) < -1e-12):
            raise SystemExit(
                f"Hit@k fell as k grew for {gene}, which is impossible for a nested "
                f"shortlist: {dict(zip(rows['k'], observed))}")
        if gene == "ALL":
            continue
        pool = int(rows["candidate_pool_size"].iloc[0])
        for _, row in rows.iterrows():
            expected = min(row["k"] / pool, 1.0)
            if abs(row["null_hit_rate"] - expected) > 1e-9:
                raise SystemExit(
                    f"{gene} null at k={int(row['k'])} is {row['null_hit_rate']:.6f}, not "
                    f"min(k/N, 1) = {expected:.6f} for a pool of {pool}")
            if row["k"] >= pool and abs(row["observed_hit_rate"] - 1.0) > 1e-9:
                raise SystemExit(
                    f"{gene} has k={int(row['k'])} >= pool {pool}, so every variant is "
                    f"inside the shortlist by construction, yet Hit@k is "
                    f"{row['observed_hit_rate']:.6f}")
    per_gene_pool = unit.groupby("gene")["candidate_pool_size"].nunique()
    if (per_gene_pool != 1).any():
        raise SystemExit(
            "the candidate pool is not constant within a gene, so a single gene-level null "
            f"is not defined:\n{unit.groupby('gene')['candidate_pool_size'].unique()}")
    n_distinct = unit.groupby("gene")["variant"].nunique()
    declared = summary.set_index("scope")["n_variants"]
    for gene, n in n_distinct.items():
        if int(declared[gene]) != int(n):
            raise SystemExit(f"{gene}: summary says {declared[gene]} distinct variants, the "
                             f"unit table holds {n}")


def main() -> None:
    args = parse_args()
    base = resolve_base(args.base)
    out_dir = reject_repo_results(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    add_harness_to_path(base)

    import harness as H  # noqa: E402  (resolved by add_harness_to_path)

    cand = {(g, s): (H.split_vars(g, s)[0] + H.split_vars(g, s)[1])
            for g in H.GENES for s in SPLITS}
    gene_cells = {g: (lambda X, lab: (X, np.asarray(lab)))(*H.load_gene(g)) for g in H.GENES}

    acc, _pools = collect(H, gene_cells, cand)
    unit = per_unit_frame(acc)
    note = verify_against_committed(unit)
    print(f"provenance check: {note}")

    columns = [f"hit_{k}" for k in KS] + [f"null_{k}" for k in KS] + ["pds", "mean_rank"]
    summary = summarize_variant(unit, columns)

    pools = unit.groupby("gene")["candidate_pool_size"].first()
    rows = []
    for _, s in summary.iterrows():
        scope = s["scope"]
        pool = int(pools[scope]) if scope != "ALL" else -1
        for k in KS:
            rows.append(dict(
                scope=scope, k=k, n_variants=int(s["n_variants"]),
                candidate_pool_size=pool,
                observed_hit_rate=float(s[f"hit_{k}"]),
                ci_lo=float(s[f"hit_{k}_lo"]),
                ci_hi=float(s[f"hit_{k}_hi"]),
                null_hit_rate=float(s[f"null_{k}"]),
                mean_rank=float(s["mean_rank"]),
                pds=float(s["pds"])))
    table = pd.DataFrame(rows)

    # Checked at full precision, then rounded for the file. The other way round, the
    # null test compares a value already rounded to six places against its closed form
    # and fails on the rounding rather than on the harness.
    check(unit, summary, table)
    ROUNDED = {"observed_hit_rate": 6, "ci_lo": 6, "ci_hi": 6, "null_hit_rate": 6,
               "mean_rank": 4, "pds": 6}
    table = table.round(ROUNDED)

    table.to_csv(out_dir / "split_half_topk.csv", index=False)
    unit.to_csv(out_dir / "split_half_topk_per_variant.csv", index=False)
    print(table.to_string(index=False))
    print(f"\nsaved -> {out_dir / 'split_half_topk.csv'}")
    print(f"saved -> {out_dir / 'split_half_topk_per_variant.csv'} ({len(unit)} rows)")


if __name__ == "__main__":
    main()
