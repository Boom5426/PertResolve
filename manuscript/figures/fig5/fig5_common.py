"""Shared style, loaders and export helpers for the refined Fig. 5 panels."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import nm_style as S  # noqa: E402
MM = 1.0 / 25.4

INK = "#202124"
SLATE = "#586773"
DARK_SLATE = "#2E3742"
MID_SLATE = "#74828D"
LIGHT_SLATE = "#A7B1B9"
GREY = "#777A7E"
MID_GREY = "#A8ADB2"
LIGHT_GREY = "#E4E7E9"
PALE_GREY = "#F3F4F5"
AMBER = "#D8A23D"

GENE_ORDER = ["TP53", "KRAS", "GATA1", "JAK1"]
GENE_COLORS = {
    "TP53": "#5185C0",
    "KRAS": "#E99D4E",
    "GATA1": "#8281B9",
    "JAK1": "#55966B",
}
HEADS = ["Ridge", "Lasso", "RF", "GBoost", "KNN", "MLP"]


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


def title_block(fig, title: str, subtitle: str | None = None, x: float = 0.02,
                *, title_y: float = 0.965, subtitle_y: float = 0.895) -> None:
    """Panel heading, wrapped to the canvas and stacked from the top.

    Both strings are wrapped rather than trusted to fit: these headings are
    prose, and prose written for a wide canvas runs off a narrow one. The
    subtitle is placed under the rendered title instead of at a fixed offset,
    so a title that wraps to two lines pushes it down rather than colliding.
    """
    title = S.wrap_to_width(fig, title, fontsize=7.4, weight="bold", x=x)
    node = fig.text(x, title_y, title, ha="left", va="top",
                    fontsize=7.4, fontweight="bold")
    if subtitle:
        renderer = fig.canvas.get_renderer()
        bottom = node.get_window_extent(renderer=renderer).y0 / fig.bbox.height
        subtitle = S.wrap_to_width(fig, subtitle, fontsize=5.7, x=x)
        fig.text(x, min(subtitle_y, bottom - 0.012), subtitle, ha="left",
                 va="top", fontsize=5.7, color=GREY)


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


def is_inhouse(method: str) -> bool:
    return any(str(method).startswith(head + "-") for head in HEADS)


def per_gene_head_means() -> pd.DataFrame:
    data = read_csv("unified_results5.csv")
    data = data[data["method"].map(is_inhouse)].copy()
    return (
        data.groupby(["gene", "method"], as_index=False)["PDS_cos"]
        .mean()
        .sort_values(["gene", "PDS_cos", "method"])
        .reset_index(drop=True)
    )


def deterministic_offsets(n: int, span: float) -> np.ndarray:
    if n <= 1:
        return np.zeros(n)
    base = np.linspace(-span, span, n)
    order = np.r_[np.arange(0, n, 2), np.arange(1, n, 2)]
    return base[order]


def mix(c0: str, c1: str, amount: float):
    a, b = to_rgb(c0), to_rgb(c1)
    return tuple(a[index] + (b[index] - a[index]) * amount for index in range(3))


def rounded_box(
    ax,
    x0: float,
    y0: float,
    width: float,
    height: float,
    *,
    edge=LIGHT_GREY,
    face="white",
    lw=0.6,
    radius=1.0,
    zorder=1,
):
    patch = FancyBboxPatch(
        (x0, y0),
        width,
        height,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        linewidth=lw,
        edgecolor=edge,
        facecolor=face,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def arrow(ax, start, end, *, color=GREY, lw=0.65, mutation=5.0, zorder=4):
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=mutation,
        linewidth=lw,
        color=color,
        shrinkA=0,
        shrinkB=0,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def mini_profile(
    ax,
    x0: float,
    x1: float,
    baseline: float,
    values,
    *,
    color=SLATE,
    amplitude=2.3,
    lw=0.85,
    zorder=3,
):
    values = np.asarray(values, dtype=float)
    scale = max(float(np.max(np.abs(values))), 1e-9)
    xs = np.linspace(x0, x1, len(values), endpoint=False)
    step = (x1 - x0) / len(values)
    ax.plot([x0, x1], [baseline, baseline], color=LIGHT_GREY, lw=0.45, zorder=zorder - 1)
    for x, value in zip(xs + step / 2.0, values):
        ax.plot(
            [x, x],
            [baseline, baseline + amplitude * value / scale],
            color=color,
            lw=lw,
            solid_capstyle="butt",
            zorder=zorder,
        )


def point_size(n: float, minimum: float = 16.0, maximum: float = 54.0) -> float:
    """Area scaling used for candidate-set size in panels f/g."""
    lo, hi = 6.0, 400.0
    value = (np.sqrt(max(n, lo)) - np.sqrt(lo)) / (np.sqrt(hi) - np.sqrt(lo))
    return float(minimum + np.clip(value, 0, 1) * (maximum - minimum))
