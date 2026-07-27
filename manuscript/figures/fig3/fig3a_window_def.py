"""Figure 3a - Definition of the detection window (schematic, drawn at final size).

Source file: none (concept panel, no data). Numbers reproduced: none; this panel
carries no measured quantity, only the definitions used by Fig. 3b-h.

One-line message: the detection window is the ratio R = D_self / D_null of the
within-variant split-half distance (replicate noise) to the variant-to-wild-type
distance (signal); R << 1 is an open window, R approximately 1 is the noise floor.

Drawn natively at its composite placement size (50 x 43 mm), so the composite scale
is ~1.0 and every glyph lands at 5.5-7 pt on the page. Deliberately monochrome
(ink / grey only): the four gene hues are reserved for the data panels, so a
neutral schematic cannot collide with a gene colour. The two-tone coding matches
panel b exactly: D_self is the emphasised quantity (ink here, gene hue in b) and
D_null is grey in both, so the same symbol is never coded two different ways.

Run:  python fig3a_window_def.py  ->  fig3a_window_def.pdf (+ .png)
Vector check:  pdfimages -list fig3a_window_def.pdf | tail -n +3 | wc -l   -> 0
"""
from __future__ import annotations

import os
import sys

import numpy as np
from matplotlib.patches import Ellipse, FancyArrowPatch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

W_MM, H_MM = 47.5, 54.5
RNG = np.random.default_rng(7)

FS_TITLE = 5.9
FS_SYM = 6.8
FS_SMALL = 5.5
FS_EQ = 7.0


def cloud(ax, cx, cy, w, h, *, wild_type=False, n=13):
    """A cell cloud: filled ellipse plus a few cells, no gene colour."""
    ec = S.GREY if wild_type else S.INK
    fc = "white" if wild_type else S.LIGHT_GREY
    ax.add_patch(Ellipse((cx, cy), w, h, facecolor=fc, edgecolor=ec, lw=0.5,
                         alpha=0.9 if wild_type else 0.55,
                         linestyle=(0, (2.2, 1.4)) if wild_type else "solid",
                         zorder=2))
    t = RNG.uniform(0, 2 * np.pi, n)
    r = np.sqrt(RNG.uniform(0, 0.62, n))
    ax.scatter(cx + r * np.cos(t) * w / 2, cy + r * np.sin(t) * h / 2,
               s=1.6, color=S.GREY if wild_type else S.INK,
               alpha=0.55, linewidths=0, zorder=3)


def dist_arrow(ax, x0, x1, y, color):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="<->",
                                 mutation_scale=4.5, lw=0.7, color=color,
                                 shrinkA=0, shrinkB=0, zorder=4))


def main() -> None:
    S.apply_rcparams()
    fig, ax = S.panel(W_MM, H_MM)
    ax.set_xlim(0, 50)
    ax.set_ylim(-0.6, 54)
    ax.axis("off")
    ax.set_position([0, 0, 1, 1])

    bw, bh = 8.4, 7.4          # cloud size, mm
    xl, xr = 8.0, 22.6         # left / right cloud centres
    x_sym = 29.2               # x where the distance symbol starts

    # ---------------- band 1: within-variant split half ----------------
    y1 = 43.5
    ax.text(0.5, 52.0, "Within variant: split into halves", fontsize=FS_TITLE,
            color=S.INK, ha="left", va="baseline", fontweight="bold")
    cloud(ax, xl, y1, bw, bh)
    cloud(ax, xr, y1, bw, bh)
    dist_arrow(ax, xl + bw / 2 + 0.5, xr - bw / 2 - 0.5, y1, S.INK)
    ax.text(xl, y1 - bh / 2 - 2.0, "half 1", fontsize=FS_SMALL, color=S.GREY,
            ha="center", va="baseline")
    ax.text(xr, y1 - bh / 2 - 2.0, "half 2", fontsize=FS_SMALL, color=S.GREY,
            ha="center", va="baseline")
    ax.text(x_sym, y1 + 1.0, r"$D_\mathrm{self}$", fontsize=FS_SYM, color=S.INK,
            ha="left", va="baseline")
    ax.text(x_sym, y1 - 3.0, "replicate noise", fontsize=FS_SMALL, color=S.GREY,
            ha="left", va="baseline")

    # ---------------- band 2: variant versus wild type ----------------
    y2 = 22.5
    ax.text(0.5, 31.0, "Variant versus wild type", fontsize=FS_TITLE,
            color=S.INK, ha="left", va="baseline", fontweight="bold")
    cloud(ax, xl, y2, bw, bh)
    cloud(ax, xr, y2, bw, bh, wild_type=True, n=9)
    dist_arrow(ax, xl + bw / 2 + 0.5, xr - bw / 2 - 0.5, y2, S.GREY)
    ax.text(xl, y2 - bh / 2 - 2.0, "variant", fontsize=FS_SMALL, color=S.GREY,
            ha="center", va="baseline")
    ax.text(xr, y2 - bh / 2 - 2.0, "wild type", fontsize=FS_SMALL, color=S.GREY,
            ha="center", va="baseline")
    ax.text(x_sym, y2 + 1.0, r"$D_\mathrm{null}$", fontsize=FS_SYM, color=S.GREY,
            ha="left", va="baseline")
    ax.text(x_sym, y2 - 3.0, "variant signal", fontsize=FS_SMALL, color=S.GREY,
            ha="left", va="baseline")

    # ---------------- band 3: the window ratio ----------------
    ax.plot([0.5, 49.5], [12.0, 12.0], color=S.LIGHT_GREY, lw=0.6, zorder=1)
    ax.text(0.5, 8.2, r"$R = D_\mathrm{self}/D_\mathrm{null}$",
            fontsize=FS_EQ, color=S.INK, ha="left", va="baseline")
    # Wording, not math symbols: Liberation Sans has no glyph for \ll or \approx,
    # so mathtext silently fell back to STIXGeneral and mixed a second font family
    # into the panel. Plain words stay in the house sans and read better at 5.5 pt.
    ax.text(0.5, 4.3, r"$R$ well below 1: detection window open", fontsize=FS_SMALL,
            color=S.INK, ha="left", va="baseline")
    ax.text(0.5, 0.8, r"$R$ near 1: at the replicate-noise floor",
            fontsize=FS_SMALL, color=S.INK, ha="left", va="baseline")

    S.save(fig, "fig3a_window_def")
    print(f"drawn at {W_MM} x {H_MM} mm (scale ~1.0 in the composite)")


if __name__ == "__main__":
    main()
