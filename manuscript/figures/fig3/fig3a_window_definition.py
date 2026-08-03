"""Fig. 3a — definition of the descriptive detection-window ratio."""
from __future__ import annotations

import numpy as np
from matplotlib.patches import Ellipse, FancyArrowPatch, FancyBboxPatch

import fig3_common as S

# Publication contract inherited from fig3_common: Arial/Helvetica sans-serif;
# svg.fonttype="none"; pdf.fonttype=42; exports .svg, .pdf, .png and .tiff at dpi=600.
# Final assembled figure target: width_mm = 183.
W_MM, H_MM = 95.0, 41.0

# The cards occupy the left of the panel and the definitions a column on the
# right. Stacked underneath, the six definition lines had 8 mm of height between
# them and three pairs overlapped. Row 1 of this figure was only 127.5 mm of the
# 183 mm page, so the panel widens rather than the type shrinking.
CARD_X = 0.560   # maps the cards' original 0-1 coordinates into the left column
DEF_X = 0.635    # left edge of the definition column


def cell_cloud(ax, center, radius, *, wild_type=False, seed_shift=0):
    cx, cy = center
    edge = S.MID_GREY if wild_type else S.INK
    fill = "white" if wild_type else S.PALE_GREY
    ax.add_patch(
        Ellipse(
            center,
            2 * radius * CARD_X,
            2 * radius,
            facecolor=fill,
            edgecolor=edge,
            lw=0.7,
            linestyle=(0, (2.2, 1.5)) if wild_type else "solid",
            zorder=2,
        )
    )
    # Deterministic schematic placement only; these dots are not observations.
    order = np.arange(12)
    theta = np.linspace(0.25, 2 * np.pi + 0.25, 12, endpoint=False)
    theta = theta + 0.43 * seed_shift
    rad = np.sqrt(0.10 + 0.50 * ((order * 7 + seed_shift) % 12) / 11.0) * radius
    ax.scatter(
        cx + rad * CARD_X * np.cos(theta),
        cy + rad * np.sin(theta),
        s=3.0,
        color=S.MID_GREY if wild_type else S.GREY,
        alpha=0.8,
        linewidths=0,
        zorder=3,
    )


def comparison_card(ax, y0, y1, title, right_label, wild_type=False):
    ax.add_patch(
        FancyBboxPatch(
            (0.035 * CARD_X, y0),
            0.93 * CARD_X,
            y1 - y0,
            boxstyle="round,pad=0.008,rounding_size=0.015",
            facecolor="white",
            edgecolor=S.LIGHT_GREY,
            lw=0.7,
            zorder=0,
        )
    )
    ax.text(
        0.065 * CARD_X,
        y1 - 0.055,
        title,
        ha="left",
        va="top",
        fontsize=6.5,
        fontweight="bold",
    )
    yc = y0 + 0.43 * (y1 - y0)
    cell_cloud(ax, (0.23 * CARD_X, yc), 0.078, seed_shift=0 if not wild_type else 2)
    cell_cloud(
        ax,
        (0.50 * CARD_X, yc),
        0.078,
        wild_type=wild_type,
        seed_shift=1 if not wild_type else 3,
    )
    ax.add_patch(
        FancyArrowPatch(
            (0.315 * CARD_X, yc),
            (0.415 * CARD_X, yc),
            arrowstyle="<->",
            mutation_scale=6,
            lw=0.8,
            color=S.INK if not wild_type else S.GREY,
            zorder=4,
        )
    )
    left_name = "half 1" if not wild_type else "variant"
    right_name = "half 2" if not wild_type else "wild type"
    ax.text(0.23 * CARD_X, yc - 0.105, left_name, ha="center", va="top", fontsize=5.7, color=S.GREY)
    ax.text(0.50 * CARD_X, yc - 0.105, right_name, ha="center", va="top", fontsize=5.7, color=S.GREY)
    ax.text(0.66 * CARD_X, yc + 0.025, right_label, ha="left", va="center", fontsize=7.2)
    ax.text(
        0.66 * CARD_X,
        yc - 0.055,
        "replicate noise" if not wild_type else "variant–WT distance",
        ha="left",
        va="center",
        fontsize=5.7,
        color=S.GREY,
    )


def main() -> None:
    S.apply_style()
    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    comparison_card(
        ax,
        0.555,
        0.960,
        "1  Within variant: split into halves",
        r"$D_\mathrm{self}$",
        wild_type=False,
    )
    comparison_card(
        ax,
        0.075,
        0.480,
        "2  Variant versus wild type",
        r"$D_\mathrm{null}$",
        wild_type=True,
    )

    ax.add_patch(
        FancyBboxPatch(
            (DEF_X - 0.022, 0.055),
            1.0 - DEF_X - 0.001,
            0.855,
            boxstyle="round,pad=0.006,rounding_size=0.012",
            facecolor=S.PALE_GREY,
            edgecolor="none",
        )
    )
    ax.text(DEF_X, 0.845, r"$R=D_\mathrm{self}/D_\mathrm{null}$",
            ha="left", va="center", fontsize=7.4, fontweight="bold")
    ax.text(DEF_X, 0.735, r"$R\approx1$", ha="left", va="center", fontsize=5.9)
    ax.text(DEF_X, 0.660, "window closed", ha="left", va="center",
            fontsize=5.9, color=S.GREY)
    ax.text(DEF_X, 0.560, r"$R$ well below 1", ha="left", va="center", fontsize=5.9)
    ax.text(DEF_X, 0.485, "window open", ha="left", va="center",
            fontsize=5.9, color=S.GREY)
    ax.text(DEF_X, 0.360, "formal detection", ha="left", va="center",
            fontsize=5.1, fontweight="bold", color=S.INK)
    ax.text(DEF_X, 0.290, r"rankability: $S>W$", ha="left", va="center",
            fontsize=5.1, fontweight="bold", color=S.INK)
    ax.text(DEF_X, 0.195, r"$S$: med $D_\mathrm{null}-$med $D_\mathrm{self}$",
            ha="left", va="center", fontsize=5.0, color=S.GREY)
    ax.text(DEF_X, 0.125, r"$W$: $Q_{97.5}(D_\mathrm{self})$", ha="left",
            va="center", fontsize=5.0, color=S.GREY)
    ax.text(DEF_X, 0.085, r"    $-\,Q_{2.5}(D_\mathrm{self})$", ha="left",
            va="center", fontsize=5.0, color=S.GREY)
    S.save(fig, "fig3a_window_definition")


if __name__ == "__main__":
    main()
