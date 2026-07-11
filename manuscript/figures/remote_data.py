"""Shared loaders for the held Fig 2/3/4/5 panels (data now local).

Dual canonical source, matching how the manuscript's numbers were produced:
  - per-METHOD forests (Fig 2b PDS, Fig 4g external models): the unified multi-seed /
    bootstrap-CI summaries in results/_remote/unified/ (PDS 0.49-0.52, Table 1 values).
  - per-GENE / per-SPLIT / per-VARIANT / per-FEATURE breakdowns (Fig 2c/d/f/g,
    Fig 4a/b/d/e/f): results/results_v4_exttheta.csv (the committed grid that reproduces
    the manuscript's breakdown numbers: gap 0.31/0.19/0.13, GATA1 Low-N PDS 0.19, 15/17).
  - Fig 3g/h, Fig 5b/c/g: the unified metrology CSVs (pairwise, classifier, oracle,
    controlled_recovery).
"""
from __future__ import annotations

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nm_style as S

ROOT = S.repo_root()
REM = ROOT / "results" / "_remote" / "unified"

# split1..6 role names (split4 = cross-gene, excluded from the main grid)
SPLIT_LABELS = {"split1": "Random", "split2": "Positional", "split3": "Mechanistic",
                "split5": "Low-depth", "split6": "Compatibility"}
SPLIT_ORDER = ["split1", "split2", "split3", "split5", "split6"]

EXTERNAL = ["scGen", "scVIDR", "Biolord", "CellFlow", "PerturbNet"]
REFS = ["Gene-mean", "WT-null"]
HEADS = ["Ridge", "Lasso", "RF", "GBoost", "KNN", "MLP"]  # x {theta, esm, esm+theta}


def feature_space(method: str) -> str:
    if method.endswith("esm+theta"):
        return "ESM+theta"
    if method.endswith("esm"):
        return "ESM"
    if method.endswith("theta"):
        return "theta"
    return "other"


def is_inhouse(method: str) -> bool:
    """The 18 feature-model heads (not external SOTA, not references)."""
    return any(method.startswith(h + "-") for h in HEADS)


# ---- breakdown grid (per variant) -----------------------------------------
def exttheta() -> pd.DataFrame:
    return pd.read_csv(ROOT / "results" / "results_v4_exttheta.csv")


# ---- per-method summaries (unified) ---------------------------------------
def definitive() -> pd.DataFrame:
    """method, PDS, ci_lo, ci_hi, crosses (bootstrap CI; multi-seed harness)."""
    return pd.read_csv(REM / "definitive_summary.csv")


def multiseed() -> pd.DataFrame:
    """method, PDS_mean, PDS_sd, lo, hi."""
    return pd.read_csv(REM / "unified_multiseed.csv")


# ---- Fig 3/5 metrology ----------------------------------------------------
def oracle() -> pd.DataFrame:
    """scope, PDS_oracle, ci_lo, ci_hi, n_var_scored (native-depth per-gene)."""
    return pd.read_csv(REM / "oracle_ceiling.csv")


def pairwise() -> pd.DataFrame:
    """gene, n_var, frac_identifiable, frac_pairs_resolvable, median_nn_dist, median_Dself, ..."""
    return pd.read_csv(REM / "pairwise_resolvability.csv")


def classifier() -> pd.DataFrame:
    """gene, detect_med_auroc, ident_med_auroc, perm_med_auroc, ..."""
    return pd.read_csv(REM / "classifier_two_sample.csv")


def controlled() -> pd.DataFrame:
    """gene, depth_m, n_var, ceiling_pds, pds_a0, pds_a1, P_correct_order, P_winner, mean_tau."""
    return pd.read_csv(REM / "controlled_recovery.csv")
