"""Fig. 2a — one held-out prediction, two non-equivalent model scores.

Direction recovery compares a prediction with the observed response of the
same held-out allele. Allele discrimination ranks that same prediction against
the full same-gene candidate set, C_g = train(g) union test(g). PDS is the
tie-aware percentile rank of the own observed allele, not a binary nearest-
neighbour indicator. Only held-out queries are scored.

Static QA contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42; outputs
.svg, .pdf, .png and .tiff at dpi=600. Final assembled figure target:
width_mm = 183.
"""
from __future__ import annotations

import os
import sys

import matplotlib as mpl
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

HERE = os.path.dirname(os.path.abspath(__file__))
W_MM, H_MM = 59.0, 63.0
PANEL_RCPARAMS = {
    "font.size": 5.0,
}

INK, GREY, LGREY = S.INK, S.GREY, S.LIGHT_GREY
TEAL = S.RESOLUTION
TEAL_FACE = S.RESOLUTION_LIGHT
CARD_FILL = "#F7F8F8"

FS_HEAD = 6.1
FS_BODY = 5.6
FS_SMALL = 5.0

PRED = np.array([0.62, -0.28, 0.86, 0.22, -0.68, 0.43, -0.16, 0.61])
MEAS = np.array([0.52, -0.38, 0.74, 0.34, -0.55, 0.33, -0.25, 0.50])


def rbox(ax, x, y, w, h, *, edge=LGREY, face="white", lw=0.6,
         radius=1.0, z=2):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        linewidth=lw, edgecolor=edge, facecolor=face, zorder=z,
    )
    ax.add_patch(patch)
    return patch


def arrow(ax, p0, p1, *, color=INK, lw=0.55, ms=4.2, rad=0.0, z=5):
    patch = FancyArrowPatch(
        p0, p1, arrowstyle="-|>", mutation_scale=ms, lw=lw, color=color,
        shrinkA=0, shrinkB=0, connectionstyle=f"arc3,rad={rad}", zorder=z,
    )
    ax.add_patch(patch)
    return patch


def stems(ax, x0, x1, yc, half_h, vals, color, *, lw=0.7, base=True):
    xs = np.linspace(x0, x1, len(vals))
    if base:
        ax.plot([x0 - 0.25, x1 + 0.25], [yc, yc], color=LGREY, lw=0.45, zorder=2)
    ax.vlines(xs, yc, yc + vals * half_h, color=color, lw=lw, zorder=3)


def main() -> None:
    S.apply_rcparams()
    mpl.rcParams.update(PANEL_RCPARAMS)
    fig, ax = S.panel(W_MM, H_MM)
    ax.set_position([0, 0, 1, 1])
    ax.set_xlim(0, W_MM)
    ax.set_ylim(0, H_MM)
    ax.set_axis_off()
    ax.add_patch(Rectangle((0, 0), W_MM, H_MM, facecolor="white",
                           edgecolor="none", zorder=-10))


    # Compact shared prediction path. Its small footprint keeps the scoring
    # distinction, rather than the already-defined workflow, visually primary.
    rbox(ax, 1.4, 52.5, 13.4, 6.0, edge=INK, lw=0.65)
    ax.text(8.1, 56.5, "held-out query", ha="center", va="center",
            fontsize=FS_SMALL, color=GREY)
    ax.text(8.1, 54.3, "TP53 R175H", ha="center", va="center", fontsize=FS_SMALL)

    arrow(ax, (15.2, 55.5), (17.0, 55.5))
    rbox(ax, 17.4, 52.5, 19.2, 6.0, edge=INK, lw=0.65)
    ax.text(27.0, 56.65, r"features  $\theta$ / ESM", ha="center", va="center",
            fontsize=FS_SMALL, color=GREY)
    ax.text(27.0, 54.85, "feature-to-response", ha="center", va="center",
            fontsize=FS_SMALL)
    ax.text(27.0, 53.55, "model", ha="center", va="center", fontsize=FS_SMALL)

    arrow(ax, (37.0, 55.5), (39.0, 55.5))
    rbox(ax, 39.4, 52.5, 18.2, 6.0, edge=INK, lw=0.65)
    ax.text(48.5, 56.8, "predicted response", ha="center", va="center",
            fontsize=FS_SMALL)
    ax.text(48.5, 55.25, r"$\hat{\delta}_v$", ha="center", va="center",
            fontsize=FS_SMALL)
    stems(ax, 43.0, 54.0, 53.75, 0.65, PRED, INK, lw=0.62)

    ax.text(29.5, 50.35, "score the same prediction", ha="center", va="center",
            fontsize=FS_SMALL, color=GREY)
    ax.plot([48.5, 48.5], [52.2, 47.7], color=GREY, lw=0.6, zorder=4)
    ax.plot([14.6, 44.4], [47.7, 47.7], color=GREY, lw=0.6, zorder=4)
    ax.plot([14.6, 14.6], [47.7, 43.0], color=GREY, lw=0.6, zorder=4)
    ax.plot([44.4, 48.5], [47.7, 47.7], color=GREY, lw=0.6, zorder=4)
    arrow(ax, (14.6, 43.0), (14.6, 42.2), color=GREY, lw=0.6)
    arrow(ax, (44.4, 47.7), (44.4, 42.2), color=GREY, lw=0.6)

    left = (1.4, 4.8, 26.4, 37.0)
    right = (31.2, 4.8, 26.4, 37.0)
    for x, y, w, h in (left, right):
        rbox(ax, x, y, w, h, edge="#B8BEC3", face="white", lw=0.7, radius=1.2)
        ax.add_patch(Rectangle((x + 0.2, y + h - 8.5), w - 0.4, 8.3,
                               facecolor=CARD_FILL, edgecolor="none", zorder=2.2))

    # Supporting direction score.
    lx, ly, lw, lh = left
    ax.text(lx + lw / 2, ly + lh - 2.0, "Direction recovery",
            ha="center", va="center", fontsize=FS_HEAD, fontweight="bold")
    ax.text(lx + lw / 2, ly + lh - 4.35, "SUPPORTING AXIS",
            ha="center", va="center", fontsize=FS_SMALL, color=GREY,
            fontweight="bold")
    ax.text(lx + lw / 2, ly + lh - 6.25,
            "prediction vs own", ha="center", va="center",
            fontsize=FS_SMALL, color=GREY)
    ax.text(lx + lw / 2, ly + lh - 7.65,
            "observed response", ha="center", va="center",
            fontsize=FS_SMALL, color=GREY)

    stems(ax, lx + 3.7, lx + 18.0, ly + 24.0, 1.7, PRED, INK, lw=0.7)
    stems(ax, lx + 3.7, lx + 18.0, ly + 18.8, 1.7, MEAS, GREY, lw=0.7)
    ax.text(lx + 19.3, ly + 24.0, "pred.", ha="left", va="center",
            fontsize=FS_SMALL, color=INK)
    ax.text(lx + 19.3, ly + 18.8, "obs.", ha="left", va="center",
            fontsize=FS_SMALL, color=GREY)

    rbox(ax, lx + 6.4, ly + 10.4, 13.6, 4.0, edge=INK, lw=0.7)
    ax.text(lx + 13.2, ly + 12.4, r"Pearson-$\delta$", ha="center",
            va="center", fontsize=FS_HEAD, fontweight="bold")
    ax.text(lx + 13.2, ly + 6.6, "Does the prediction recover\nthe observed direction?",
            ha="center", va="center", fontsize=FS_SMALL, linespacing=1.2)

    # Primary allele-ranking score. The illustrative order deliberately places
    # the own allele second and uses no distance arrows, so no success is implied.
    rx, ry, rw, rh = right
    ax.text(rx + rw / 2, ry + rh - 2.0, "Allele discrimination",
            ha="center", va="center", fontsize=FS_HEAD, fontweight="bold")
    ax.text(rx + rw / 2, ry + rh - 4.35, "PRIMARY MODEL SCORE",
            ha="center", va="center", fontsize=FS_SMALL, color=TEAL,
            fontweight="bold")
    ax.text(rx + rw / 2, ry + rh - 6.15, "same-gene candidates",
            ha="center", va="center", fontsize=FS_SMALL, color=GREY)
    ax.text(rx + rw / 2, ry + rh - 7.65,
            # "union" spelled out, not \cup: U+222A is absent from every
            # Arial-metric sans, so mathtext silently takes that one glyph from
            # STIXGeneral and puts a second typeface in the panel.
            r"$C_g$ = train($g$) union test($g$)", ha="center", va="center",
            fontsize=FS_SMALL, color=GREY)

    ax.text(rx + rw / 2, ry + 27.0, r"ranked by similarity to $\hat{\delta}_v$",
            ha="center", va="center", fontsize=FS_SMALL, color=GREY)
    rank_rows = [
        (ry + 23.7, "1", "R248W", "sibling", False),
        (ry + 19.8, "2", "R175H", "own observed", True),
        (ry + 15.9, "3", "R273C", "sibling", False),
    ]
    for cy, rank, allele, role, own in rank_rows:
        edge = TEAL if own else LGREY
        face = TEAL_FACE if own else "white"
        ax.text(rx + 2.1, cy, rank, ha="center", va="center", fontsize=FS_SMALL,
                color=TEAL if own else GREY, fontweight="bold")
        rbox(ax, rx + 3.5, cy - 1.45, 21.2, 2.9, edge=edge, face=face,
             lw=0.7 if own else 0.5, radius=0.7)
        if own:
            ax.text(rx + 14.1, cy + 0.42, allele, ha="center", va="center",
                    fontsize=FS_SMALL, color=TEAL, fontweight="bold")
            ax.text(rx + 14.1, cy - 0.72, role, ha="center", va="center",
                    fontsize=FS_SMALL, color=TEAL)
        else:
            ax.text(rx + 7.0, cy, allele, ha="left", va="center",
                    fontsize=FS_SMALL, color=INK)
            ax.text(rx + 23.5, cy, role, ha="right", va="center",
                    fontsize=FS_SMALL, color=GREY)

    rbox(ax, rx + 3.6, ry + 8.4, 19.2, 5.0, edge=TEAL, face="white", lw=0.75)
    ax.text(rx + 13.2, ry + 10.9, "PDS", ha="center", va="center",
            fontsize=FS_HEAD, color=TEAL, fontweight="bold")
    ax.text(rx + 13.2, ry + 5.6, "Tie-aware percentile rank of\nthe own observed allele",
            ha="center", va="center", fontsize=FS_SMALL, linespacing=1.2)
    ax.text(rx + 13.2, ry + 2.0, "chance = 0.50", ha="center", va="center",
            fontsize=FS_SMALL, color=TEAL, fontweight="bold")

    ax.text(W_MM / 2, 2.2,
            "Held-out queries only; candidate set contains all same-gene variants.",
            ha="center", va="center", fontsize=FS_SMALL, color=GREY)

    S.save(fig, os.path.join(HERE, "fig2a_evaluation"), exact=True,
           formats=("pdf", "png", "svg", "tiff"))


if __name__ == "__main__":
    main()
