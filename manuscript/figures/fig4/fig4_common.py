"""Shared style, loaders and export helpers for the refined Fig. 4 panels."""
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

INK = "#202124"
SLATE = "#586773"
PDS_BLUE = "#5B7FA3"
DIRECTION_BLUE = "#3F6F9F"
NEGATIVE = "#B77A67"
AMBER = "#D8A23D"
GREY = "#777A7E"
MID_GREY = "#A8ADB2"
LIGHT_GREY = "#E4E7E9"
PALE_GREY = "#F3F4F5"

GENE_ORDER = ["TP53", "KRAS", "GATA1", "JAK1"]
GENE_COLORS = {
    "TP53": "#5185C0",
    "KRAS": "#E99D4E",
    "GATA1": "#8281B9",
    "JAK1": "#55966B",
}

HEADS = ["Ridge", "Lasso", "RF", "GBoost", "KNN", "MLP"]
FEATURE_ORDER = ["theta", "esm", "esm+theta"]
FEATURE_LABELS = {
    "theta": r"$\theta$",
    "esm": "ESM",
    "esm+theta": r"ESM+$\theta$",
}
SPLIT_ORDER = ["split1", "split2", "split3", "split5", "split6"]
SPLIT_LABELS = {
    "split1": "Random",
    "split2": "Positional",
    "split3": "Mechanistic",
    "split5": "Low-depth",
    "split6": "Compatibility",
}
SPLIT_GENE_COUNTS = {
    "split1": 4,
    "split2": 4,
    "split3": 4,
    "split5": 2,
    "split6": 3,
}


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


def breakdown() -> pd.DataFrame:
    return read_csv("results_v4_exttheta.csv")


def is_inhouse(method: str) -> bool:
    return any(str(method).startswith(head + "-") for head in HEADS)


def feature_space(method: str) -> str:
    if method.endswith("esm+theta"):
        return "esm+theta"
    if method.endswith("esm"):
        return "esm"
    if method.endswith("theta"):
        return "theta"
    return "other"


def head_name(method: str) -> str:
    return method.split("-", 1)[0]


def ordered_methods() -> list[str]:
    return [f"{head}-{feature}" for feature in FEATURE_ORDER for head in HEADS]


def wilson_interval(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n <= 0:
        return math.nan, math.nan
    p = k / n
    denominator = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denominator
    half = z * math.sqrt(
        p * (1.0 - p) / n + z * z / (4.0 * n * n)
    ) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def gene_ticklabels(ax, axis: str = "x") -> None:
    ticks = ax.get_xticklabels() if axis == "x" else ax.get_yticklabels()
    for tick, gene in zip(ticks, GENE_ORDER):
        tick.set_color(GENE_COLORS[gene])
        tick.set_fontweight("bold")


def deterministic_offsets(n: int, span: float) -> np.ndarray:
    if n <= 1:
        return np.zeros(n)
    return np.linspace(-span, span, n)
