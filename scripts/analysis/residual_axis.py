#!/usr/bin/env python3
"""Score every model on a second axis: discrimination after the gene-shared programme.

Most of a variant's measured response is the programme its gene shares with every sibling.
A prediction that reproduces only that programme is right about the gene and silent about
the allele, and the Discussion of the manuscript already says so: predictors "should be
judged on whether they recover allele-specific residual structure, not the gene-shared
response that dominates direction-based scores". Nothing in the paper scores that way. This
does.

**The residual axis is added, never substituted.** Both axes are computed for every method
here, from the same predictions, the same held-out variants and the same candidate sets. No
model is retrained and no existing definition changes.

**Every score is checked against a permutation null computed on the same axis.** Removing a
large component shared by every candidate changes the geometry a cosine sees, and a new axis
on which several methods appear above chance, in a study whose result is that they are not,
has to be assumed to be an artefact until it survives a null built the same way. Predictions
are permuted across the held-out variants of a gene, which destroys the correspondence
between a prediction and its target while preserving every marginal property of both, and
the whole scoring is repeated. A score is reported as informative only where it exceeds that
null, not merely where it exceeds 0.5.

**The gene mean comes from training variants only.** The earlier decomposition in
``results/reviewer_controls/residual_decomp_*.csv`` centred on a mean that included the
variants being scored: on the measurement side over every variant, so each variant's
residual lost a 1/n share of its own signal and the residuals acquired a spurious mutual
anti-correlation that a ranking statistic benefits from; on the model side over the test
variants alone, so the definition of the residual saw the held-out set. Neither is fatal at
these sample sizes, but neither is defensible in a table that a reader will treat as a model
score, so the mean here is computed over each split's training variants and nothing else.

Emits one row per (method, gene, split) and a pooled summary, with both axes side by side.

Usage:
  python residual_axis.py --out OUT_DIR [--base BASE]
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from alleleperturb.paths import (
    add_harness_to_path,
    reject_repo_results,
    require_inputs,
    resolve_base,
)

_ap = argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--base", default=None,
                 help="VCCompass compute workspace holding unified/ (env: VCCOMPASS_BASE)")
_ap.add_argument("--out", required=True, metavar="OUT_DIR",
                 help="directory receiving residual_axis_models.csv and "
                      "residual_axis_summary.csv; may not be inside the repository's results/")
_ap.add_argument("--n-boot", type=int, default=2000,
                 help="bootstrap resamples over held-out variants (default: %(default)s)")
_ap.add_argument("--n-perm", type=int, default=50,
                 help="permutations of predictions across held-out variants, per method "
                      "(default: %(default)s)")
_args = _ap.parse_args()

BASE = resolve_base(_args.base)
OUT_DIR = reject_repo_results(_args.out)
OUT_DIR.mkdir(parents=True, exist_ok=True)
require_inputs(BASE / "unified" / "real_deltas.npz", BASE / "allele_perturb_bench.csv")
add_harness_to_path(BASE)

import harness as H  # noqa: E402

SPLITS = ["split1", "split2", "split3", "split5", "split6"]
NBOOT = _args.n_boot
N_PERM = _args.n_perm
np.random.seed(0)

REAL = {k: v for k, v in np.load(BASE / "unified" / "real_deltas.npz").items()}


def tie_aware_pds(pred: np.ndarray, target_index: int, truths: np.ndarray) -> float:
    """Rank one prediction against every candidate's measured profile, ties at mid-rank."""
    if np.linalg.norm(pred) < 1e-12:
        return 0.5
    pn = pred / np.linalg.norm(pred)
    tn = truths / (np.linalg.norm(truths, axis=1, keepdims=True) + 1e-12)
    dist = 1.0 - tn @ pn
    target = dist[target_index]
    n = len(dist)
    if n < 2:
        return 0.5
    less = int((dist < target - 1e-12).sum())
    eq = int((np.abs(dist - target) <= 1e-12).sum())
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1)


def score_split(gene: str, split: str, predictions, *, permute_seed: int | None = None) -> list[dict]:
    """Both axes for every held-out variant of one gene and split.

    The candidate set is the split's training and held-out variants together, matching the
    convention every other score in this study uses. The gene mean removed for the residual
    axis is computed over the training variants only, so the held-out variants take no part
    in defining the quantity they are scored on.
    """
    train, test = H.split_vars(gene, split)
    candidates = [v for v in train + test if f"{gene}__{v}" in REAL]
    held_out = [v for v in test if v in candidates
                and f"{gene}__{split}__{v}" in predictions.files]
    train_present = [v for v in train if f"{gene}__{v}" in REAL]
    if len(candidates) < 3 or not held_out or len(train_present) < 2:
        return []

    truths = np.stack([REAL[f"{gene}__{v}"] for v in candidates]).astype(np.float64)
    gene_mean = np.stack([REAL[f"{gene}__{v}"] for v in train_present]).mean(axis=0)
    residual_truths = truths - gene_mean
    index = {v: i for i, v in enumerate(candidates)}

    # Under the null, each held-out variant is scored against another held-out variant's
    # prediction. Every prediction and every truth is the same as before; only which goes
    # with which changes, so anything the axis produces here comes from its geometry rather
    # than from a model knowing something.
    source = dict(zip(held_out, held_out))
    if permute_seed is not None and len(held_out) > 1:
        shuffled = list(held_out)
        rng = np.random.RandomState(permute_seed)
        for _ in range(20):
            rng.shuffle(shuffled)
            if all(a != b for a, b in zip(held_out, shuffled)):
                break
        source = dict(zip(held_out, shuffled))

    rows = []
    for v in held_out:
        pred = np.asarray(predictions[f"{gene}__{split}__{source[v]}"], dtype=np.float64)
        rows.append(dict(
            gene=gene, split=split, variant=v,
            pds=tie_aware_pds(pred, index[v], truths),
            residual_pds=tie_aware_pds(pred - gene_mean, index[v], residual_truths),
            n_candidates=len(candidates), n_train_in_mean=len(train_present),
        ))
    return rows


def oracle_rows(gene: str, split: str) -> list[dict]:
    """The same two axes for a second measurement of each variant, not for a model.

    This is the ceiling both axes are read against. It uses the identical candidate set and
    the identical training-only gene mean, so the comparison between a model and the ceiling
    is like for like on each axis.
    """
    train, test = H.split_vars(gene, split)
    candidates = [v for v in train + test if f"{gene}__{v}" in REAL]
    train_present = [v for v in train if f"{gene}__{v}" in REAL]
    if len(candidates) < 3 or len(train_present) < 2:
        return []

    X, lab = H.load_gene(gene)
    lab = np.asarray(lab)
    wt_tags = ("WT", "wt", "WT_control")
    rows = []
    for seed in range(15):
        rng = np.random.RandomState(1000 + 7 * seed + H._GENE_SEED[gene])
        wt = np.where(np.isin(lab, wt_tags))[0]
        rng.shuffle(wt)
        half = max(len(wt) // 2, 1)
        wm_a, wm_b = X[wt[:half][:300]].mean(0), X[wt[half:][:300]].mean(0)
        first, second = {}, {}
        for v in candidates:
            idx = np.where(lab == v)[0]
            if len(idx) < 5:
                continue
            rng.shuffle(idx)
            h = len(idx) // 2
            first[v] = X[idx[:h][:300]].mean(0) - wm_a
            second[v] = X[idx[h:][:300]].mean(0) - wm_b
        usable = [v for v in candidates if v in first]
        train_here = [v for v in train_present if v in first]
        held_out = [v for v in test if v in first]
        if len(usable) < 3 or len(train_here) < 2 or not held_out:
            continue
        truths = np.stack([second[v] for v in usable]).astype(np.float64)
        gene_mean = np.stack([second[v] for v in train_here]).mean(axis=0)
        residual_truths = truths - gene_mean
        index = {v: i for i, v in enumerate(usable)}
        for v in held_out:
            query = np.asarray(first[v], dtype=np.float64)
            rows.append(dict(
                gene=gene, split=split, variant=v, seed=seed,
                pds=tie_aware_pds(query, index[v], truths),
                residual_pds=tie_aware_pds(query - gene_mean, index[v], residual_truths),
            ))
    return rows


def bootstrap_mean(values: np.ndarray, seed: int = 0) -> tuple[float, float, float]:
    """Mean of a per-variant score with a percentile interval over variants."""
    values = np.asarray(values, dtype=np.float64)
    rng = np.random.RandomState(seed)
    n = len(values)
    draws = np.array([values[rng.randint(0, n, n)].mean() for _ in range(NBOOT)])
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return float(values.mean()), float(lo), float(hi)


print("=== the replicate ceiling on both axes ===")
ceiling_rows = []
for gene in H.GENES:
    per_variant = {}
    for split in SPLITS:
        for row in oracle_rows(gene, split):
            key = (row["variant"], row["split"])
            per_variant.setdefault(key, []).append(row)
    if not per_variant:
        continue
    full = np.array([np.mean([r["pds"] for r in rs]) for rs in per_variant.values()])
    resid = np.array([np.mean([r["residual_pds"] for r in rs]) for rs in per_variant.values()])
    f, f_lo, f_hi = bootstrap_mean(full)
    r, r_lo, r_hi = bootstrap_mean(resid)
    ceiling_rows.append(dict(gene=gene, n_scored=len(full),
                             ceiling_pds=round(f, 3), ceiling_lo=round(f_lo, 3),
                             ceiling_hi=round(f_hi, 3),
                             ceiling_residual_pds=round(r, 3),
                             ceiling_residual_lo=round(r_lo, 3),
                             ceiling_residual_hi=round(r_hi, 3)))
    print(f"  {gene:6s} n={len(full):4d}  PDS {f:.3f} [{f_lo:.3f}, {f_hi:.3f}]   "
          f"residual-PDS {r:.3f} [{r_lo:.3f}, {r_hi:.3f}]")
ceiling = pd.DataFrame(ceiling_rows)
ceiling.to_csv(OUT_DIR / "residual_axis_ceiling.csv", index=False)

print("\n=== every method on both axes ===")
methods = sorted(os.path.basename(f)[:-4]
                 for f in glob.glob(str(BASE / "unified" / "preds5" / "*.npz")))
per_variant_rows, summary_rows, gene_rows = [], [], []
for method in methods:
    predictions = np.load(BASE / "unified" / "preds5" / f"{method}.npz")
    rows = [r for gene in H.GENES for split in SPLITS
            for r in score_split(gene, split, predictions)]
    if not rows:
        print(f"  {method:20s} no scorable held-out variant")
        continue
    for row in rows:
        row["method"] = method
    per_variant_rows.extend(rows)

    frame = pd.DataFrame(rows)
    # average the repeated scorings of a variant across splits first, so a variant held out
    # by several splits does not count several times in the interval
    per_var = frame.groupby(["gene", "variant"])[["pds", "residual_pds"]].mean()
    f, f_lo, f_hi = bootstrap_mean(per_var["pds"].to_numpy())
    r, r_lo, r_hi = bootstrap_mean(per_var["residual_pds"].to_numpy())

    null_full, null_resid = [], []
    for perm in range(N_PERM):
        prows = [x for gene in H.GENES for split in SPLITS
                 for x in score_split(gene, split, predictions, permute_seed=1000 + perm)]
        if not prows:
            continue
        pv = pd.DataFrame(prows).groupby(["gene", "variant"])[["pds", "residual_pds"]].mean()
        null_full.append(pv["pds"].mean())
        null_resid.append(pv["residual_pds"].mean())
    p_full = (np.sum(np.array(null_full) >= f) + 1) / (len(null_full) + 1) if null_full else np.nan
    p_resid = ((np.sum(np.array(null_resid) >= r) + 1) / (len(null_resid) + 1)
               if null_resid else np.nan)

    summary_rows.append(dict(method=method, n_variants=len(per_var),
                             pds=round(f, 3), pds_lo=round(f_lo, 3), pds_hi=round(f_hi, 3),
                             residual_pds=round(r, 3), residual_lo=round(r_lo, 3),
                             residual_hi=round(r_hi, 3),
                             null_pds=round(float(np.mean(null_full)), 3) if null_full else np.nan,
                             null_residual_pds=round(float(np.mean(null_resid)), 3) if null_resid else np.nan,
                             perm_p_pds=round(float(p_full), 4),
                             perm_p_residual=round(float(p_resid), 4),
                             pds_crosses_chance=bool(f_lo <= 0.5 <= f_hi),
                             residual_crosses_chance=bool(r_lo <= 0.5 <= r_hi)))
    print(f"  {method:20s} n={len(per_var):4d}  PDS {f:.3f} (null {np.mean(null_full):.3f}, "
          f"P={p_full:.3f})   residual {r:.3f} (null {np.mean(null_resid):.3f}, P={p_resid:.3f})")

    for gene, sub in frame.groupby("gene"):
        pv = sub.groupby("variant")[["pds", "residual_pds"]].mean()
        gene_rows.append(dict(method=method, gene=gene, n_variants=len(pv),
                              pds=round(float(pv["pds"].mean()), 3),
                              residual_pds=round(float(pv["residual_pds"].mean()), 3)))

pd.DataFrame(per_variant_rows).to_csv(OUT_DIR / "residual_axis_per_variant.csv", index=False)
summary = pd.DataFrame(summary_rows).sort_values("pds", ascending=False)
summary.to_csv(OUT_DIR / "residual_axis_models.csv", index=False)

print(f"\nsaved -> {OUT_DIR}")
pd.DataFrame(gene_rows).to_csv(OUT_DIR / "residual_axis_per_gene.csv", index=False)
print(f"\nmethods above their own permutation null on the residual axis (P < 0.05): "
      f"{int((summary.perm_p_residual < 0.05).sum())} of {len(summary)}")
print(f"\nmethods whose PDS interval excludes chance: "
      f"{int((~summary.pds_crosses_chance).sum())} of {len(summary)}")
print(f"methods whose residual-PDS interval excludes chance: "
      f"{int((~summary.residual_crosses_chance).sum())} of {len(summary)}")
