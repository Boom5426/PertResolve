"""Figure 3a - Split-half definition of the measurement window (schematic).

Programmatic (vector) concept panel, no data. A variant's cells are split into two
random halves; their distance is the within-variant noise D_self, the distance to
wild-type is the signal D_null, and R = D_self/D_null defines the window. Two
mini-scenes contrast an open window (WT far -> large D_null -> R << 1) with a closed
window (WT near -> D_null approximately D_self -> R approximately 1). Kept
programmatic (not AI) so it stays vector and palette-consistent.

Run:  python fig3a_schematic.py  ->  fig3a_schematic.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

from matplotlib.patches import Ellipse, FancyArrowPatch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

VAR = S.GENE_COLORS["TP53"]      # generic variant colour (blue)
VAR_EC = "#3E6A9E"
WT = "#B9B9B9"
NOISE = S.GENE_COLORS["KRAS"]    # orange for D_self


def blob(ax, x, y, w, h, fc, ec, label=None, fs=5.5, tc="white"):
    ax.add_patch(Ellipse((x, y), w, h, facecolor=fc, edgecolor=ec, lw=0.6, alpha=0.9))
    if label:
        ax.text(x, y, label, ha="center", va="center", fontsize=fs, color=tc)


def arrow(ax, p0, p1, color):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="<->", mutation_scale=6,
                                 lw=0.9, color=color, shrinkA=1, shrinkB=1))


def main() -> None:
    S.apply_rcparams()
    fig, ax = S.panel(52, 46)
    ax.set_xlim(0, 10); ax.set_ylim(0, 9.4); ax.axis("off")

    # ---------- top: construction ----------
    blob(ax, 1.6, 7.4, 2.0, 1.5, WT, "#8F8F8F", "WT", tc="#444")
    blob(ax, 6.0, 8.0, 1.7, 1.15, VAR, VAR_EC, "half 1")
    blob(ax, 6.7, 6.7, 1.7, 1.15, VAR, VAR_EC, "half 2")

    arrow(ax, (2.6, 7.5), (5.3, 7.6), VAR)          # D_null
    ax.text(3.9, 8.35, r"$D_\mathrm{null}$", ha="center", fontsize=7, color=VAR,
            fontweight="bold")
    ax.text(3.9, 6.95, "signal (variant vs WT)", ha="center", fontsize=4.6, color=S.GREY)

    arrow(ax, (6.15, 7.45), (6.55, 7.25), NOISE)    # D_self
    ax.text(7.9, 7.4, r"$D_\mathrm{self}$", ha="left", fontsize=7, color=NOISE,
            fontweight="bold")
    ax.text(7.9, 6.9, "noise (half vs half)", ha="left", fontsize=4.6, color=S.GREY)

    ax.text(5.0, 5.35, r"window ratio   $R = D_\mathrm{self}\,/\,D_\mathrm{null}$",
            ha="center", va="center", fontsize=7, color=S.INK)

    # ---------- divider ----------
    ax.plot([0.3, 9.7], [4.55, 4.55], color="#E0E0E0", lw=0.6)

    # ---------- bottom-left: open window (WT far -> big D_null) ----------
    blob(ax, 0.9, 3.0, 1.0, 0.8, WT, "#8F8F8F")
    blob(ax, 3.5, 3.1, 0.9, 0.7, VAR, VAR_EC)
    blob(ax, 3.8, 2.7, 0.9, 0.7, VAR, VAR_EC)
    arrow(ax, (1.45, 3.0), (2.95, 3.0), VAR)
    ax.text(2.3, 1.9, r"open window: $R \ll 1$", ha="center", fontsize=5.6,
            color=S.GENE_COLORS["JAK1"], fontweight="bold")
    ax.text(2.3, 1.35, "measurable perturbation", ha="center", fontsize=4.6, color=S.GREY)

    # ---------- bottom-right: closed window (WT near -> D_null ~ D_self) ----------
    blob(ax, 6.2, 3.0, 1.0, 0.8, WT, "#8F8F8F")
    blob(ax, 7.4, 3.1, 0.9, 0.7, VAR, VAR_EC)
    blob(ax, 7.7, 2.7, 0.9, 0.7, VAR, VAR_EC)
    arrow(ax, (6.75, 3.0), (7.0, 3.0), VAR)
    ax.text(7.4, 1.9, r"closed window: $R \approx 1$", ha="center", fontsize=5.6,
            color=S.GENE_COLORS["GATA1"], fontweight="bold")
    ax.text(7.4, 1.35, "little usable window", ha="center", fontsize=4.6, color=S.GREY)

    S.save(fig, "fig3a_schematic")
    print("drawn schematic")


if __name__ == "__main__":
    main()
