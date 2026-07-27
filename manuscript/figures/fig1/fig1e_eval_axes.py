"""Figure 1e - AllelePerturb-Eval separates three evaluation axes.

Replaces the raster panel ``Fig1e_Eval_highres.pdf`` with a vector panel drawn
at its FINAL composite size (76 x 35 mm, placed at 76 mm, scale 1.0), so all
type lands at 5.5-6.5 pt on the page.

Source file: none (metric-definition panel). No score, benchmark value or
result appears here; only metric names taken verbatim from the Fig. 1e caption
and the Results text: direction recovery (Pearson-δ, δ-cosine), allele
discrimination (PDS) and differential-expression fidelity (DE overlap, DE-LFC
rank correlation, direction agreement).

One-line message: evaluation is split into direction recovery, allele
discrimination and DE fidelity, and allele discrimination is the defining axis.

Run:  python fig1e_eval_axes.py
Out:  fig1e_eval_axes.pdf (+ .png preview)
Vector check:  pdfimages -list fig1e_eval_axes.pdf | tail -n +3 | wc -l  -> 0
"""
from __future__ import annotations

import os
import sys

import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

W, H = 76.0, 35.0  # final placement size in mm

COLS = [
    dict(cx=13.0, header="Direction recovery",
         metrics=["Pearson-δ", "δ-cosine"],
         sub="does the prediction\npoint the right way?"),
    dict(cx=38.0, header="Allele discrimination",
         metrics=["PDS", "d(own) < d(others)"],
         sub="does it identify its\nown allele?"),
    dict(cx=63.0, header="DE fidelity",
         metrics=["DE overlap", "DE-LFC rank corr.", "direction agreement"],
         sub="do the changed\ngenes agree?"),
]

Y_HEADER = 33.6
Y_METRIC = [18.6, 16.3, 14.0]
Y_SUB = 12.0
BANNER_H = 6.2


def icon_direction(ax, cx):
    """Two roughly parallel profiles: prediction and observation point the same way."""
    x = np.linspace(cx - 10.5, cx + 3.0, 60)
    t = np.linspace(0, 1, 60)
    base = 1.9 * np.sin(2 * np.pi * t * 1.3) + 1.1 * np.sin(2 * np.pi * t * 0.55)
    y0 = 26.3
    y_solid = y0 + base
    y_dash = y0 - 3.2 + base * 0.85
    ax.plot(x, y_solid, color=S.INK, lw=0.7)
    ax.plot(x, y_dash, color=S.GREY, lw=0.7, ls=(0, (2, 1.4)))
    ax.text(cx + 3.8, y_solid[-1], "pred.", ha="left", va="center", fontsize=5.5,
            color=S.INK)
    ax.text(cx + 3.8, y_dash[-1], "obs.", ha="left", va="center", fontsize=5.5,
            color=S.GREY)


def icon_discrimination(ax, cx):
    """One prediction, ranked against sibling variant profiles; own variant is nearest."""
    ax.add_patch(FancyBboxPatch((cx - 9.4, 23.0), 4.4, 3.4,
                                boxstyle="round,pad=0,rounding_size=0.7",
                                facecolor="white", edgecolor=S.INK, lw=0.5))
    ax.text(cx - 7.2, 24.7, "δ̂", ha="center", va="center", fontsize=6.0,
            color=S.INK)
    ax.add_patch(FancyArrowPatch((cx - 4.6, 24.7), (cx - 1.8, 24.7),
                                 arrowstyle="-|>", mutation_scale=4.5, lw=0.5,
                                 color=S.GREY, shrinkA=0, shrinkB=0))
    widths = [7.6, 5.4, 6.2, 4.6]
    for i, w in enumerate(widths):
        y = 28.4 - i * 2.6
        hot = i == 0
        ax.add_patch(Rectangle((cx - 1.2, y - 0.85), w,
                               1.7, facecolor=S.HOTSPOT if hot else S.LIGHT_GREY,
                               edgecolor="none"))
    ax.text(cx + 7.2, 28.4, "own", ha="left", va="center", fontsize=5.5,
            color=S.HOTSPOT)
    ax.text(cx + 5.0, 21.0, "siblings", ha="left", va="center", fontsize=5.5,
            color=S.GREY)


def icon_de(ax, cx):
    """Predicted vs measured changed-gene lists, one gene shared between them."""
    for j, (x0, col, lab) in enumerate([(cx - 8.6, S.INK, "pred."),
                                        (cx + 1.4, S.GREY, "obs.")]):
        ax.text(x0 + 3.5, 29.6, lab, ha="center", va="bottom", fontsize=5.5,
                color=col)
        for i in range(4):
            y = 27.8 - i * 2.2
            shared = (j == 0 and i in (0, 2)) or (j == 1 and i in (1, 2))
            ax.add_patch(Rectangle((x0, y - 0.7), 7.0, 1.4,
                                   facecolor=S.HOTSPOT if shared else S.LIGHT_GREY,
                                   edgecolor="none"))


ICONS = [icon_direction, icon_discrimination, icon_de]


def main() -> None:
    S.apply_rcparams()
    fig, ax = S.panel(W, H)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(Rectangle((0, 0), W, H, facecolor="white", edgecolor="none",
                           zorder=0))

    # thin separators between the three axes
    for x in (25.5, 50.5):
        ax.plot([x, x], [BANNER_H + 1.4, 34.6], color=S.LIGHT_GREY, lw=0.5)

    for col, draw_icon in zip(COLS, ICONS):
        cx = col["cx"]
        ax.text(cx, Y_HEADER, col["header"], ha="center", va="center",
                fontsize=6.0, fontweight="bold", color=S.INK)
        draw_icon(ax, cx)
        for y, m in zip(Y_METRIC, col["metrics"]):
            ax.text(cx, y, m, ha="center", va="center", fontsize=5.8, color=S.INK)
        ax.text(cx, Y_SUB, col["sub"], ha="center", va="top", fontsize=5.5,
                color=S.GREY, linespacing=1.35)

    # ---- bottom banner: the central question -------------------------------
    # a hairline rule, not a filled box: the tint used before was off-palette
    # and read as a decorative panel box
    ax.plot([0.6, W - 0.6], [BANNER_H + 0.1] * 2, color=S.LIGHT_GREY, lw=0.5,
            solid_capstyle="butt")
    ax.text(W / 2, BANNER_H / 2 - 0.3,
            "Central question: can models distinguish alleles, not merely\n"
            "recover the shared direction of a perturbed gene?",
            ha="center", va="center", fontsize=5.5, color=S.INK, linespacing=1.35)

    S.save(fig, "fig1e_eval_axes")


if __name__ == "__main__":
    main()
