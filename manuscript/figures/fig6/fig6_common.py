"""Shared visual language, data paths and export helpers for refined Fig. 6."""
# Final assembled figure target: width_mm = 183; independent panel canvases are
# intentionally exported at their working sizes before journal-layout assembly.
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import pandas as pd

import nm_style as S  # noqa: E402

# Source tables live in the repository, not in a private copy beside the panels.
# A panel that reads its own copy keeps plotting the old numbers after the
# committed table is regenerated.
RAW = S.repo_root() / "results"
DERIVED = S.repo_root() / "results" / "fig6_derived"
MM = 1.0 / 25.4

INK = "#202124"
DARK_SLATE = "#2E3742"
MID_SLATE = "#5F6B76"
SLATE = "#75828D"
LIGHT_SLATE = "#9AA7B3"
GREY = "#777A7E"
MID_GREY = "#A8ADB2"
LIGHT_GREY = "#DDE1E4"
PALE_GREY = "#F3F4F5"
BLUE_PALE = "#E8F0F7"
AMBER = "#D8A23D"
AMBER_PALE = "#F7EEDC"

DATASET_COLORS = {
    "TP53": "#5185C0",
    "KRAS": "#E99D4E",
    "GATA1": "#8281B9",
    "JAK1": "#55966B",
    "Replogle": DARK_SLATE,
    "Norman": MID_SLATE,
    "Adamson": LIGHT_SLATE,
}
GENE_ORDER = ["TP53", "KRAS", "GATA1", "JAK1"]
EXTERNAL_ORDER = ["Replogle", "Norman", "Adamson"]


def read_csv(name: str):
    """Read a committed canonical table by name, as the other figures do.

    Panels in this figure historically read only from ``DERIVED``. Panels that draw from
    ``results/canonical/`` need the same resolver the rest of the figures use, which
    searches the committed data directories and raises rather than falling back.
    """
    return pd.read_csv(S.data_path(name))


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


def figure(width_mm: float, height_mm: float):
    return plt.figure(figsize=(width_mm * MM, height_mm * MM))


def panel(width_mm: float, height_mm: float):
    return plt.subplots(figsize=(width_mm * MM, height_mm * MM))


def title_block(
    fig,
    title: str,
    subtitle: str | None = None,
    *,
    x: float = 0.02,
    title_y: float = 0.965,
    subtitle_y: float = 0.895,
) -> None:
    """Panel heading, wrapped to the canvas and stacked from the top.

    Both strings are wrapped rather than trusted to fit: these headings are
    prose, and prose written for a wide canvas runs off a narrow one. The
    subtitle is placed under the rendered title instead of at a fixed offset,
    so a title that wraps to two lines pushes it down rather than colliding.
    """
    title = S.wrap_to_width(fig, title, fontsize=7.5, weight="bold", x=x)
    node = fig.text(x, title_y, title, ha="left", va="top",
                    fontsize=7.5, fontweight="bold")
    if subtitle:
        renderer = fig.canvas.get_renderer()
        bottom = node.get_window_extent(renderer=renderer).y0 / fig.bbox.height
        subtitle = S.wrap_to_width(fig, subtitle, fontsize=5.7, x=x)
        fig.text(x, min(subtitle_y, bottom - 0.012), subtitle, ha="left",
                 va="top", fontsize=5.7, color=GREY)


def vector_gradient(ax, cmap, extent, *, steps: int = 96, zorder: float = 2) -> None:
    """Re-exported from nm_style so panels reach it through this module."""
    S.vector_gradient(ax, cmap, extent, steps=steps, zorder=zorder)


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


def rounded_box(
    ax,
    x: float,
    y: float,
    width: float,
    height: float,
    *,
    edge=LIGHT_GREY,
    face="white",
    lw=0.65,
    radius=1.1,
    zorder=1,
):
    patch = FancyBboxPatch(
        (x, y),
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


def arrow(ax, start, end, *, color=GREY, lw=0.65, mutation=5.5, zorder=4):
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
