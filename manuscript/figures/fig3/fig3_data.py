"""Committed-data loaders for Figure 3, faithfully replicating the canonical
native-depth rankability aggregation from scripts/figures/fig_config.py so the
panels here are self-contained (no dependency on the legacy figure code/style).

Canonical definition (locked): native max-depth (deepest split-half bin per
perturbation) + all perturbations + direct split-half energy distance, PCA-50.
Verified 2026-07-11 to reproduce canonical_numbers.json exactly:
  D_self/D_null mean  TP53 0.965 / KRAS 1.004 / GATA1 0.878 / JAK1 0.210
  rankable (native)   TP53 0%    / KRAS 0%    / GATA1 2.4%  / JAK1 90%
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

RESULTS = S.repo_root() / "results"


def _rank_table() -> pd.DataFrame:
    return pd.read_csv(RESULTS / "second_probe_rankability_table.csv")


def native_rankability(dataset: str, metric: str = "edist", space: str = "pca") -> pd.DataFrame:
    """Per-perturbation rows at each perturbation's deepest split-half bin."""
    df = _rank_table()
    sub = df[(df.dataset == dataset) & (df.metric == metric) & (df.space == space)]
    if len(sub) == 0:
        return sub
    native = sub.loc[sub.groupby("perturbation")["n_work"].idxmax()].copy()
    native["ratio"] = native["D_self"] / native["D_null"].replace(0, np.nan)
    return native


def rankable_fraction(dataset: str) -> tuple[float, int]:
    """Return (rankable_fraction_0to1, n) at native depth."""
    n = native_rankability(dataset)
    return n["rankable"].mean(), len(n)


def power_curve(space: str = "pca50") -> pd.DataFrame:
    df = pd.read_csv(RESULTS / "split_half_power_curve.csv")
    return df[df.space == space].copy()


def dself_dnull_ci() -> dict:
    """Bootstrap CIs for the per-gene D_self/D_null aggregate (Fig 3 quoted values)."""
    return json.load(open(RESULTS / "bootstrap_CIs.json"))["dself_dnull_ci"]


def unrankable_canonical() -> dict:
    return json.load(open(RESULTS / "unrankable_canonical.json"))["canonical_native"]
