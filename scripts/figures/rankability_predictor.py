#!/usr/bin/env python
"""Honest per-perturbation rankability predictor (leave-one-dataset-out).

This replaces the earlier ``second_probe_predictor_results.csv`` (headline
"AUROC 0.974"), which a 5-way adversarial audit showed to be a
config-identity + pseudo-replication artifact: it was trained ROW-level over the
~18 correlated (space x metric x n_work) rows per perturbation, and its
discrimination came from evaluation-configuration identity (which metric/space a
row used) rather than from any pilot-estimable perturbation feature. Under a
genuine per-perturbation setup, TP53 and KRAS are uniformly unrankable, so an
AUROC for them is mathematically undefined.

This script predicts per-perturbation rankability from pilot-estimable,
non-leaky features under leave-one-dataset-out (LODO) cross-validation.

Label + aggregation follow the paper's canonical rankability definition (one row
per perturbation at its deepest available split-half bin; energy distance;
PCA-50 space; matches ``fig_config.native_rankability``).

Excluded by construction, with reasons:
  * S, W, ratio, D_self, D_null   -- ``rankable == (S > W)`` exactly, so these
    ARE the label (target leakage);
  * space, metric, level, n_work  -- evaluation-configuration identity, not a
    property of the perturbation, and the driver of the previous artifact;
  * "PCA explained-variance ratio" and "gene-space dimensionality" -- named in an
    earlier Methods draft but absent from the data (the latter is a per-dataset
    constant, inert for within-fold ranking anyway).

Held-out datasets whose label is single-class (all-unrankable) are reported as
NOT EVALUABLE (AUROC undefined); they are never silently dropped or imputed.

Usage:
    python scripts/figures/rankability_predictor.py \
        --table results/second_probe_rankability_table.csv \
        --out   results/rankability_predictor_honest.csv
"""
from __future__ import annotations

import argparse
import json
from typing import Sequence

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

LEAKY_COLUMNS = ("S", "W", "ratio", "D_self", "D_null")  # rankable == (S > W)
CONFIG_IDENTITY = ("space", "metric", "level", "n_work")  # not per-perturbation
SEED = 0
N_BOOT = 2000
DATASETS = ("Replogle", "Norman", "Adamson", "TP53", "KRAS", "GATA1", "JAK1")
ALLELEPERTURB = ("TP53", "KRAS", "GATA1", "JAK1")


def per_perturbation(df: pd.DataFrame, metric: str = "edist", space: str = "pca") -> pd.DataFrame:
    """Collapse to one row per (dataset, perturbation) at its deepest split-half bin.

    Matches the canonical ``fig_config.native_rankability`` definition
    (``groupby('perturbation')['n_work'].idxmax()``) used for Fig 3d/5d.
    """
    sub = df[(df["metric"] == metric) & (df["space"] == space) & (df["n_cells"] > 0)].copy()
    idx = sub.groupby(["dataset", "perturbation"])["n_work"].idxmax()
    return sub.loc[idx].reset_index(drop=True)


def _auroc_with_ci(y: np.ndarray, p: np.ndarray, seed: int = SEED, n_boot: int = N_BOOT):
    """AUROC plus percentile bootstrap 95% CI over perturbations.

    Returns (auroc, lo, hi) or (nan, nan, nan) if the label is single-class.
    """
    if len(np.unique(y)) < 2:
        return np.nan, np.nan, np.nan
    auroc = roc_auc_score(y, p)
    rng = np.random.RandomState(seed)
    boots = []
    n = len(y)
    for _ in range(n_boot):
        idx = rng.randint(0, n, n)
        if len(np.unique(y[idx])) < 2:
            continue
        boots.append(roc_auc_score(y[idx], p[idx]))
    if not boots:
        return auroc, np.nan, np.nan
    return auroc, float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def lodo(pp: pd.DataFrame, feats: Sequence[str], seed: int = SEED) -> pd.DataFrame:
    """Leave-one-dataset-out logistic regression on per-perturbation features."""
    feats = list(feats)
    rows = []
    for ho in DATASETS:
        te = pp[pp["dataset"] == ho]
        tr = pp[pp["dataset"] != ho]
        n_pert = len(te)
        rate = float(te["rankable"].mean()) if n_pert else np.nan
        evaluable = te["rankable"].nunique() > 1
        if evaluable and tr["rankable"].nunique() > 1:
            scaler = StandardScaler().fit(tr[feats])
            model = LogisticRegression(max_iter=2000, random_state=seed).fit(
                scaler.transform(tr[feats]), tr["rankable"].astype(int)
            )
            p = model.predict_proba(scaler.transform(te[feats]))[:, 1]
            auroc, lo, hi = _auroc_with_ci(te["rankable"].astype(int).values, p, seed=seed)
        else:
            auroc, lo, hi = np.nan, np.nan, np.nan
        rows.append(
            dict(
                held_out=ho,
                n_perturbations=int(n_pert),
                rankable_rate=round(rate, 4) if n_pert else np.nan,
                evaluable=bool(evaluable),
                auroc=round(auroc, 4) if not np.isnan(auroc) else np.nan,
                auroc_lo=round(lo, 4) if not np.isnan(lo) else np.nan,
                auroc_hi=round(hi, 4) if not np.isnan(hi) else np.nan,
                features="+".join(feats),
            )
        )
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", default="results/second_probe_rankability_table.csv",
                    help="canonical rankability table (input)")
    ap.add_argument("--out", default="results/rankability_predictor_honest.csv",
                    help="honest per-dataset predictor results (output)")
    ap.add_argument("--metric", default="edist")
    ap.add_argument("--space", default="pca")
    args = ap.parse_args()

    df = pd.read_csv(args.table)

    # Fail loudly if the label is not exactly (S > W): the whole audit rests on it.
    agree = (df["rankable"].astype(int) == (df["S"] > df["W"]).astype(int)).mean()
    if abs(agree - 1.0) > 1e-9:
        raise SystemExit(f"label check failed: rankable == (S>W) holds for only {agree:.4f} of rows")

    pp = per_perturbation(df, metric=args.metric, space=args.space)
    pp["log_effect"] = np.log10(pp["effect_size"].clip(lower=1e-6))
    pp["log_ncells"] = np.log10(pp["n_cells"].clip(lower=1))

    print(f"Per-perturbation table: {len(pp)} perturbations "
          f"({args.metric}/{args.space}, deepest split-half bin).")
    print("Per-dataset rankable rate (single-class => NOT EVALUABLE):")
    for ds in DATASETS:
        d = pp[pp["dataset"] == ds]
        deg = " <- degenerate (all one class)" if d["rankable"].nunique() < 2 else ""
        print(f"  {ds:9s} n={len(d):5d}  rankable={d['rankable'].mean():.3f}{deg}")

    feature_sets = {
        "effect_size": ["log_effect"],
        "effect_size+log_ncells": ["log_effect", "log_ncells"],
    }

    all_out = []
    summary = {}
    for name, feats in feature_sets.items():
        res = lodo(pp, feats)
        res["feature_set"] = name
        all_out.append(res)
        ev = res[res["evaluable"]]
        mean_all = float(ev["auroc"].mean())
        ap4 = res[res["held_out"].isin(ALLELEPERTURB) & res["evaluable"]]
        summary[name] = dict(
            evaluable_datasets=list(ev["held_out"]),
            not_evaluable=list(res.loc[~res["evaluable"], "held_out"]),
            mean_auroc_evaluable=round(mean_all, 4),
            per_dataset={r.held_out: (None if np.isnan(r.auroc) else r.auroc) for r in res.itertuples()},
            alleleperturb_evaluable=list(ap4["held_out"]),
            mean_auroc_alleleperturb_evaluable=round(float(ap4["auroc"].mean()), 4) if len(ap4) else None,
        )
        print(f"\n=== feature set: {name} ===")
        print(res[["held_out", "n_perturbations", "rankable_rate", "evaluable",
                   "auroc", "auroc_lo", "auroc_hi"]].to_string(index=False))
        print(f"  mean AUROC over EVALUABLE datasets ({len(ev)}): {mean_all:.3f}")

    out = pd.concat(all_out, ignore_index=True)
    out.to_csv(args.out, index=False)
    print(f"\nWrote {args.out}")
    print("\nSUMMARY (honest):")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
