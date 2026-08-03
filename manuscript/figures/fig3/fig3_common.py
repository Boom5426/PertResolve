"""Shared style, data loaders and export helpers for the refined Fig. 3 panels."""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import nm_style as S  # noqa: E402
MM = 1.0 / 25.4

GENE_ORDER = ["TP53", "KRAS", "GATA1", "JAK1"]
GENE_COLORS = {
    "TP53": "#5185C0",
    "KRAS": "#E99D4E",
    "GATA1": "#8281B9",
    "JAK1": "#55966B",
}
INK = "#202124"
GREY = "#777A7E"
MID_GREY = "#A8ADB2"
LIGHT_GREY = "#E7E9EB"
PALE_GREY = "#F4F5F6"


def apply_style() -> None:
    """House style from nm_style, then this figure's local overrides.

    The font stack deliberately is not restated here. A local stack that lists
    Arial first and DejaVu Sans last renders in DejaVu Sans, silently and
    without error, on any machine without Arial; and the Arial shipped by
    Debian/Ubuntu lacks U+0302, so mathtext takes every \\hat accent from
    STIXGeneral and puts a second typeface in the panel. nm_style.resolve_sans()
    picks a family by glyph coverage and raises when none qualifies.
    """
    S.apply_rcparams()
    plt.rcParams.update(
        {
            "font.size": 7.0,
            "axes.titlesize": 7.0,
            "axes.labelsize": 7.0,
            "xtick.labelsize": 6.0,
            "ytick.labelsize": 6.0,
            "legend.fontsize": 5.6,
            "axes.linewidth": 0.55,
            "lines.linewidth": 0.85,
            "patch.linewidth": 0.55,
            "xtick.major.width": 0.55,
            "ytick.major.width": 0.55,
            "text.color": INK,
            "axes.labelcolor": INK,
            "axes.edgecolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
        }
    )


def figure(w_mm: float, h_mm: float):
    return plt.figure(figsize=(w_mm * MM, h_mm * MM))


def panel(w_mm: float, h_mm: float):
    return plt.subplots(figsize=(w_mm * MM, h_mm * MM))


def despine(ax, keep=("left", "bottom")) -> None:
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(side in keep)


def save(fig, stem: str) -> None:
    """Write the declared canvas as PDF, PNG, SVG and 600 dpi LZW TIFF.

    exact=True keeps the physical canvas the panel declared; the compositor
    places every panel at the width it was drawn at, so a tight bounding box
    would silently change the on-page type size.
    """
    S.save(fig, HERE / stem, exact=True, formats=("pdf", "png", "svg", "tiff"))


def read_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(S.data_path(name))


def read_json(name: str):
    with S.data_path(name).open() as stream:
        return json.load(stream)


def rank_table() -> pd.DataFrame:
    return read_csv("second_probe_rankability_table.csv")


def native_rankability(gene: str) -> pd.DataFrame:
    df = rank_table()
    sub = df[
        (df["dataset"] == gene)
        & (df["metric"] == "edist")
        & (df["space"] == "pca")
    ]
    native = sub.loc[sub.groupby("perturbation")["n_work"].idxmax()].copy()
    native["window_ratio"] = native["D_self"] / native["D_null"].replace(0, np.nan)
    native["variant"] = native["perturbation"].str.split("__").str[-1]
    return native


def ratio_ci() -> dict:
    return read_json("bootstrap_CIs.json")["dself_dnull_ci"]


def nominal_depth_join() -> pd.DataFrame:
    """Native-depth metrology rows joined to the uncapped nominal cell counts."""
    bench = read_csv("allele_perturb_bench.csv")[["gene", "variant", "n_cells"]]
    rows = []
    for gene in GENE_ORDER:
        native = native_rankability(gene)
        merged = native.merge(
            bench,
            left_on=["dataset", "variant"],
            right_on=["gene", "variant"],
            how="left",
            suffixes=("", "_nominal"),
            validate="one_to_one",
        )
        if merged["n_cells_nominal"].isna().any():
            missing = merged.loc[merged["n_cells_nominal"].isna(), "perturbation"].tolist()
            raise ValueError(f"Nominal depth missing for {gene}: {missing}")
        rows.append(merged)
    return pd.concat(rows, ignore_index=True)


def wilson_interval(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    """Two-sided Wilson score interval for a binomial proportion."""
    if n <= 0:
        return math.nan, math.nan
    p = k / n
    den = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / den
    half = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / den
    return max(0.0, center - half), min(1.0, center + half)


def gene_ticklabels(ax, axis: str = "x") -> None:
    ticks = ax.get_xticklabels() if axis == "x" else ax.get_yticklabels()
    for tick, gene in zip(ticks, GENE_ORDER):
        tick.set_color(GENE_COLORS[gene])
        tick.set_fontweight("bold")


def size_from_nominal_cells(values) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    lo, hi = np.log10(50.0), np.log10(10000.0)
    scaled = (np.log10(np.clip(values, 50.0, 10000.0)) - lo) / (hi - lo)
    return 6.0 + 28.0 * np.clip(scaled, 0.0, 1.0)
