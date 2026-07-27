"""Committed-data loaders + shared dataset palette for Figure 6, so all panels are
self-contained and mutually consistent (no dependency on the legacy figure code).

Verified sources (2026-07-11):
  6b/6c  results/rankability_predictor_honest.csv (feature_set='effect_size')
  6d     results/pilot_validation/pilot_validation_summary.json (+ README for SNR/null)
  6e     results/second_probe_rankability_table.csv via native_rankability
  6f     results/canonical_numbers.json (unrankable native vs matched-50)
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

# dataset palette, house system only. Allele genes keep their unique house gene
# hues (TP53 blue, KRAS orange, GATA1 purple, JAK1 green) so colour means "gene"
# wherever a gene is plotted; the external gene-level atlases are given the house
# neutral slate ramp (nm_style.FEATURE_COLORS) instead of invented hues, dark to
# light in the order in which they rank across panels c/d/f. The previous
# teal/mauve/tan set collided with the gene hues (tan vs KRAS orange in
# particular) and made the figure speak two colour languages.
DATASET_COLORS = {
    "TP53": S.GENE_COLORS["TP53"], "KRAS": S.GENE_COLORS["KRAS"],
    "GATA1": S.GENE_COLORS["GATA1"], "JAK1": S.GENE_COLORS["JAK1"],
    "Replogle": S.FEATURE_COLORS["ESM+theta"],   # dark slate  #2E3742
    "Norman": S.FEATURE_COLORS["ESM"],           # mid slate   #5F6B76
    "Adamson": S.FEATURE_COLORS["theta"],        # light slate #9AA7B3
    "VCC": S.GREY, "sciPlex": S.LIGHT_GREY,
}
# datasets with non-degenerate LODO labels, in a sensible display order
LODO_DATASETS = ["Replogle", "Norman", "Adamson", "GATA1", "JAK1"]


def native_rankability(dataset: str, metric: str = "edist", space: str = "pca") -> pd.DataFrame:
    """Per-perturbation rows at each perturbation's deepest split-half bin
    (matches fig_config.native_rankability; adds a `ratio` column)."""
    df = pd.read_csv(RESULTS / "second_probe_rankability_table.csv")
    sub = df[(df.dataset == dataset) & (df.metric == metric) & (df.space == space)]
    if len(sub) == 0:
        return sub
    n = sub.loc[sub.groupby("perturbation")["n_work"].idxmax()].copy()
    n["ratio"] = n["D_self"] / n["D_null"].replace(0, np.nan)
    return n


def honest_auroc(feature_set: str = "effect_size") -> pd.DataFrame:
    df = pd.read_csv(RESULTS / "rankability_predictor_honest.csv")
    return df[df.feature_set == feature_set].reset_index(drop=True)


def pilot_summary() -> dict:
    return json.load(open(RESULTS / "pilot_validation" / "pilot_validation_summary.json"))


def canonical() -> dict:
    return json.load(open(RESULTS / "canonical_numbers.json"))
