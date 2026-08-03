"""Figure 1d (left) — auditable allele-derived representations.

The six bars are the actual corrected-table values for TP53 R175H, not decorative
bullets.  A separate ESM branch makes clear that a six-dimensional biophysical
vector and a 1,280-dimensional protein embedding are alternative inputs.  Only
the theta representation is projected in the adjacent PCA panel.

Run: python fig1d_theta_schematic.py
Output: fig1d_theta_schematic.svg/.pdf/.png

Imported export contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42;
outputs .svg, .pdf, .png and .tiff at dpi=600. Final assembled figure target:
width_mm = 183.
"""
from __future__ import annotations

import os
import sys

import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

from pathlib import Path

HERE = Path(__file__).resolve().parent

W_MM, H_MM = 69.0, 42.0
FEATURES = [
    ("Δ hydrophobicity", "d_hydro"),
    ("Δ side-chain volume", "d_vol"),
    ("Δ charge", "d_charge"),
    ("fold-core location", "fold_core"),
    ("functional switch", "cat_switch"),
    ("external hotspot", "is_hotspot"),
]


def arrow(ax, x0, y0, x1, y1, *, color=S.GREY):
    ax.add_patch(
        FancyArrowPatch(
            (x0, y0),
            (x1, y1),
            arrowstyle="-|>",
            mutation_scale=4.8,
            lw=0.55,
            color=color,
            shrinkA=0,
            shrinkB=0,
        )
    )


def main() -> None:
    S.apply_rcparams()
    data = S.load_bench(exclude_wt=True, corrected=True)
    row = data[(data.gene == "TP53") & (data.variant == "R175H")]
    if len(row) != 1:
        raise ValueError(f"Expected one corrected-table TP53 R175H row, found {len(row)}")
    values = row.iloc[0]

    fig, ax = S.panel(W_MM, H_MM)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_xlim(0, W_MM)
    ax.set_ylim(0, H_MM)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(Rectangle((0, 0), W_MM, H_MM, facecolor="white", edgecolor="none"))

    # Variant condition -----------------------------------------------------
    blue = S.GENE_COLORS["TP53"]
    ax.text(1.0, 34.8, "TP53", fontsize=6.4, fontweight="bold", color=S.INK, ha="left")
    ax.add_patch(Rectangle((1.0, 29.9), 12.5, 2.4, facecolor=S.LIGHT_GREY, edgecolor="none"))
    ax.add_patch(Rectangle((3.6, 29.9), 6.6, 2.4, facecolor=blue, alpha=0.52, edgecolor="none"))
    ax.plot([7.2, 7.2], [29.4, 33.1], color=S.HOTSPOT, lw=0.85)
    ax.scatter([7.2], [33.1], s=8, color=S.HOTSPOT, edgecolor=S.INK, linewidths=0.35)
    ax.text(7.2, 26.7, "R175H", fontsize=6.0, ha="center")
    ax.text(7.2, 23.6, "R → H", fontsize=5.4, color=S.GREY, ha="center")

    arrow(ax, 14.6, 29.8, 17.7, 29.8)
    arrow(ax, 14.6, 25.1, 17.7, 10.5)

    # Theta vector ---------------------------------------------------------
    ax.text(18.5, 39.2, r"$\theta$  6D biophysical feature vector",
            fontsize=6.1, fontweight="bold", ha="left")
    # Not "standardized values": only d_hydro, d_vol and d_charge are z-scored
    # over the 470 variants (mean 0, sd 1). fold_core, cat_switch and is_hotspot
    # are on their native 0-1 annotation scale, so the hotspot bar reads +1.00
    # where its z-score would be +2.35. Labelling the axis as the values the
    # model actually receives is both true and the more informative statement.
    ax.text(66.6, 36.8, r"$\theta$ as supplied to the model", fontsize=5.0,
            color=S.GREY, ha="right")

    x_label = 18.5
    x_zero = 48.5
    max_width = 9.0
    ax.plot([x_zero, x_zero], [15.7, 36.4], color=S.LIGHT_GREY, lw=0.55)
    ys = np.linspace(34.4, 17.2, len(FEATURES))
    for y, (label, column) in zip(ys, FEATURES):
        value = float(values[column])
        clipped = float(np.clip(value, -1.0, 1.0))
        colour = S.HOTSPOT if column == "is_hotspot" else S.FEATURE_COLORS["theta"]
        ax.text(x_label, y, label, fontsize=5.1, ha="left", va="center",
                color=S.INK)
        x1 = x_zero + clipped * max_width
        left = min(x_zero, x1)
        width = max(0.22, abs(x1 - x_zero))
        ax.add_patch(
            Rectangle(
                (left, y - 0.72),
                width,
                1.44,
                facecolor=colour,
                edgecolor="none",
                alpha=0.92 if column == "is_hotspot" else 0.75,
            )
        )
        ax.text(59.1, y, f"{value:+.2f}", fontsize=5.0,
                ha="left", va="center", color=S.GREY)

    ax.text(63.0, 14.2, "shared space", fontsize=5.0,
            color=S.RESOLUTION, ha="center")
    arrow(ax, 60.8, 12.9, 67.4, 12.9, color=S.RESOLUTION)

    # Alternative ESM representation --------------------------------------
    ax.add_patch(
        FancyBboxPatch(
            (18.5, 4.6),
            41.0,
            8.0,
            boxstyle="round,pad=0.25,rounding_size=0.8",
            facecolor="white",
            edgecolor=S.LIGHT_GREY,
            lw=0.6,
        )
    )
    ax.text(20.0, 10.5, "ESM", fontsize=5.8, fontweight="bold",
            color=S.FEATURE_COLORS["ESM"], ha="left")
    ax.text(27.2, 10.5, "mutant protein sequence", fontsize=5.0,
            color=S.GREY, ha="left")
    ax.add_patch(Rectangle((20.0, 6.4), 21.5, 2.0,
                           facecolor=S.FEATURE_COLORS["ESM"], alpha=0.18,
                           edgecolor=S.FEATURE_COLORS["ESM"], lw=0.4))
    ax.plot([30.4, 30.4], [6.1, 8.7], color=S.HOTSPOT, lw=0.85)
    arrow(ax, 42.7, 7.4, 47.5, 7.4)
    for i, alpha in enumerate((0.25, 0.42, 0.60, 0.78, 0.95)):
        ax.add_patch(Rectangle((48.3 + i * 1.8, 6.3), 1.25, 2.2,
                               facecolor=S.FEATURE_COLORS["ESM"], alpha=alpha,
                               edgecolor="none"))
    ax.text(57.7, 7.4, "1,280D", fontsize=5.0, color=S.GREY,
            ha="left", va="center")

    ax.text(1.0, 1.15, "Input representations are derived from the allele, not its response cells.",
            fontsize=5.0, color=S.GREY, ha="left")

    S.save(fig, HERE / "fig1d_theta_schematic", exact=True,
           formats=("pdf", "png", "svg", "tiff"))


if __name__ == "__main__":
    main()
