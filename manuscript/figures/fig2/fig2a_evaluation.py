"""Fig 2a: unified evaluation of one held-out variant (orientation schematic).

Source file: none. This panel plots NO benchmark numbers; it is a conceptual
schematic. The profile glyphs are fixed illustrative shapes (hard-coded arrays
below), not data, and no PDS / Pearson value, chance line or model name appears.
Content follows fig2a_prompt.md and the Fig 2 caption: a held-out variant's
predicted response is scored two ways, direction recovery (Pearson-delta) and
allele discrimination (PDS, whether it identifies its own variant among the
other held-out variants).

Numbers reproduced: none (schematic only).
One-line message: one prediction, two distinct questions.
Run:  python fig2a_evaluation.py  ->  fig2a_evaluation.pdf (+ .png)
Vector check:  pdfimages -list fig2a_evaluation.pdf | tail -n +3 | wc -l  ->  0

Drawn at its FINAL placement size (59 x 62 mm) so the composite scale is ~1.0
and every label lands at >= 5.2 pt on the page.
"""
from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

HERE = os.path.dirname(os.path.abspath(__file__))

W_MM, H_MM = 59.0, 63.0

# Colour discipline: the four gene hues are reserved for genes (Fig 2f/2g and
# Figs 1/3/6) and the slate ramp for feature spaces (Fig 2b/2c/2d). The two
# evaluation branches are therefore NOT colour-coded; they are separated by
# position and label, drawn in ink, and the single accent is the house HOTSPOT
# highlight marking the one candidate that is the variant's own profile.
INK, GREY, LGREY = S.INK, S.GREY, S.LIGHT_GREY
BRANCH = INK                      # both branches: structure, not category
HILITE = S.HOTSPOT                # highlight of the correct ("own") candidate

FS_HEAD = 6.0
FS_BODY = 5.6
FS_SMALL = 5.4

# fixed illustrative profile shapes (schematic glyphs, not measured data)
PRED = np.array([0.55, -0.30, 0.85, 0.20, -0.70, 0.40, -0.15, 0.65, -0.45, 0.25])
MEAS = np.array([0.45, -0.42, 0.70, 0.35, -0.58, 0.30, -0.28, 0.52, -0.35, 0.38])
CAND = {
    "A": np.array([0.50, -0.35, 0.78, 0.28, -0.64, 0.36, -0.20, 0.58, -0.40, 0.31]),
    "B": np.array([-0.35, 0.60, -0.20, 0.72, 0.15, -0.55, 0.40, -0.25, 0.50, -0.60]),
    "C": np.array([0.20, 0.45, -0.65, -0.15, 0.55, 0.25, -0.50, -0.30, 0.35, 0.60]),
}


def stems(ax, x0, x1, yc, half_h, vals, color, lw=0.7, base=True):
    """Draw a mini stem profile between x0 and x1, centred on yc."""
    xs = np.linspace(x0, x1, len(vals))
    if base:
        ax.plot([x0 - 0.6, x1 + 0.6], [yc, yc], color=LGREY, lw=0.5, zorder=2)
    for x, v in zip(xs, vals):
        ax.plot([x, x], [yc, yc + v * half_h], color=color, lw=lw,
                solid_capstyle="round", zorder=3)


def rbox(ax, xc, yc, w, h, edge, face="white", lw=0.6, z=2):
    b = FancyBboxPatch((xc - w / 2, yc - h / 2), w, h,
                       boxstyle="round,pad=0,rounding_size=1.0",
                       linewidth=lw, edgecolor=edge, facecolor=face, zorder=z)
    ax.add_patch(b)
    return b


def arrow(ax, p0, p1, color=INK, lw=0.5, rad=0.0, z=5, ms=4.0):
    a = FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms,
                        lw=lw, color=color, shrinkA=0, shrinkB=0,
                        connectionstyle=f"arc3,rad={rad}", zorder=z)
    ax.add_patch(a)
    return a


def main() -> None:
    S.apply_rcparams()
    # Keep Greek/maths glyphs in the same Arial-metric sans as the body text.
    # nm_style sets the text font but not mathtext, whose default (DejaVu Sans)
    # would embed a second typeface for every $\theta$, $\delta$ and subscript.
    plt.rcParams.update({"mathtext.fontset": "custom",
                         "mathtext.rm": "Liberation Sans",
                         "mathtext.it": "Liberation Sans:italic",
                         "mathtext.bf": "Liberation Sans:bold",
                         "mathtext.default": "it"})
    fig, ax = S.panel(W_MM, H_MM)
    ax.set_position([0, 0, 1, 1])
    ax.set_xlim(0, W_MM)
    ax.set_ylim(-0.5, H_MM - 0.5)
    ax.set_axis_off()
    # pin the tight bbox to the full canvas so the placed scale stays ~1.0
    ax.add_patch(Rectangle((0, -0.5), W_MM, H_MM, facecolor="white",
                           edgecolor="none", zorder=0))

    xc = 29.5

    # ---- 1. the held-out variant -----------------------------------------
    rbox(ax, xc, 57.5, 32, 8.0, edge=INK)
    ax.text(xc, 59.6, "held-out variant", ha="center", va="center",
            fontsize=FS_SMALL, color=GREY, zorder=6)
    ax.text(xc, 55.9, "TP53 R175H", ha="center", va="center",
            fontsize=FS_HEAD, color=INK, zorder=6)

    # ---- 2. features -> model --------------------------------------------
    arrow(ax, (xc, 53.5), (xc, 50.7))
    ax.text(xc + 1.6, 52.1, r"features: $\theta$ or ESM", ha="left", va="center",
            fontsize=FS_SMALL, color=GREY, zorder=6)

    rbox(ax, xc, 47.6, 40, 6.0, edge=INK)
    ax.text(xc, 47.6, "feature-to-response model", ha="center", va="center",
            fontsize=FS_BODY, color=INK, zorder=6)

    # ---- 3. predicted response -------------------------------------------
    arrow(ax, (xc, 44.6), (xc, 42.0))
    ax.text(xc, 41.0, r"predicted response  $\hat{\delta}_v$", ha="center",
            va="center", fontsize=FS_BODY, color=INK, zorder=6)
    stems(ax, xc - 12, xc + 12, 36.6, 2.9, PRED, INK, lw=0.7)

    # ---- 4. two-way branch ------------------------------------------------
    arrow(ax, (xc - 4.0, 33.4), (16.0, 32.3), color=BRANCH, lw=0.6, rad=0.15)
    arrow(ax, (xc + 4.0, 33.4), (43.0, 32.3), color=BRANCH, lw=0.6, rad=-0.15)

    # ================= LEFT: direction recovery ============================
    lx = 14.8
    ax.text(lx, 30.7, "Direction recovery", ha="center", va="center",
            fontsize=FS_HEAD, color=INK, fontweight="bold", zorder=6)
    ax.text(lx, 28.1, "predicted vs measured", ha="center", va="center",
            fontsize=FS_SMALL, color=GREY, zorder=6)

    stems(ax, 3.6, 20.4, 23.0, 2.6, PRED, INK, lw=0.7)
    stems(ax, 3.6, 20.4, 17.6, 2.6, MEAS, GREY, lw=0.7)
    ax.text(21.8, 23.0, "pred.", ha="left", va="center", fontsize=FS_SMALL,
            color=INK, zorder=6)
    ax.text(21.8, 17.6, "meas.", ha="left", va="center", fontsize=FS_SMALL,
            color=GREY, zorder=6)

    rbox(ax, lx, 11.9, 20, 4.4, edge=INK, face="white", lw=0.8)
    ax.text(lx, 11.9, r"Pearson-$\delta$", ha="center", va="center",
            fontsize=FS_HEAD, color=INK, fontweight="bold", zorder=6)
    ax.text(lx, 7.5, "Does the response\npoint the right way?", ha="center",
            va="center", fontsize=FS_SMALL, color=INK, linespacing=1.35, zorder=6)

    # ================= RIGHT: allele discrimination ========================
    rx = 43.8
    ax.text(rx, 30.7, "Allele discrimination", ha="center", va="center",
            fontsize=FS_HEAD, color=INK, fontweight="bold", zorder=6)
    ax.text(rx, 28.1, "vs other held-out variants", ha="center", va="center",
            fontsize=FS_SMALL, color=GREY, zorder=6)

    # predicted profile on the left of this column
    stems(ax, 31.5, 36.0, 21.0, 2.4, PRED, INK, lw=0.6)
    ax.text(33.7, 24.4, r"$\hat{\delta}_v$", ha="center", va="center",
            fontsize=FS_BODY, color=INK, zorder=6)

    cand_y = [25.0, 21.0, 17.0]
    cand_lab = ["R175H (own)", "R248W", "R273C"]
    cand_key = ["A", "B", "C"]
    for y, lab, key in zip(cand_y, cand_lab, cand_key):
        own = "own" in lab
        col = HILITE if own else GREY
        rbox(ax, 49.0, y, 17.4, 3.4, edge=col if own else LGREY, lw=0.6)
        stems(ax, 41.5, 45.2, y, 1.2, CAND[key], col, lw=0.5, base=False)
        ax.text(51.7, y, lab, ha="center", va="center", fontsize=FS_SMALL,
                color=col, zorder=6)
        arrow(ax, (37.0, 21.0), (40.4, y), color=col if own else LGREY,
              lw=0.6 if own else 0.45, rad=0.0, ms=3.4)

    rbox(ax, rx, 11.9, 20, 4.4, edge=INK, face="white", lw=0.8)
    ax.text(rx, 11.9, "PDS", ha="center", va="center",
            fontsize=FS_HEAD, color=INK, fontweight="bold", zorder=6)
    ax.text(rx, 7.5, "Is its own variant\nthe nearest match?", ha="center",
            va="center", fontsize=FS_SMALL, color=INK, linespacing=1.35, zorder=6)

    # ---- footer -----------------------------------------------------------
    ax.plot([4.0, 55.0], [3.4, 3.4], color=LGREY, lw=0.5, zorder=1)
    ax.text(xc, 1.8, "one prediction, two distinct questions", ha="center",
            va="center", fontsize=FS_SMALL, color=GREY, style="italic", zorder=6)

    S.save(fig, os.path.join(HERE, "fig2a_evaluation"))


if __name__ == "__main__":
    main()
