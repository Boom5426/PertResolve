"""Figure 5e - Why the measurement window governs benchmark validity (schematic).

Source file: none. Illustrative concept panel, no measured quantity is plotted and no
real PDS value, dataset or model name appears. It takes the five synthetic predictors
of panel d, whose true quality order P1 < ... < P5 is fixed by construction, and shows
the same true order (i) scrambled by a benchmark whose ceiling sits near chance and
(ii) recovered by a benchmark whose ceiling clears the floor. The quantitative version
of this statement is panels f and g.
Numbers reproduced: none (deliberate; see fig5e_prompt.md "Do NOT include").
Drawn at 85 x 36 mm, i.e. its exact placement size in fig5_assemble.tex (scale 1.0).

Run:  python fig5e_leaderboard.py  ->  fig5e_leaderboard.pdf (+ .png)
Vector check:  pdfimages -list fig5e_leaderboard.pdf | tail -n +3 | wc -l   # -> 0
"""
from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

W_MM, H_MM = 85.0, 36.0
# Same sequential slate ramp as panel d: P1 (worst, light) -> P5 (best, dark).
# Colour is redundant with the P1..P5 text, so the panel also reads in greyscale,
# and the scrambled board is visible as a broken lightness order, not only as
# crossing lines. No gene hue and no red is used: red/green are reserved for signed
# quantities, and blue/orange/purple/green mean TP53/KRAS/GATA1/JAK1 in b, c, f, g.
from matplotlib.colors import to_rgb


def _mix(c0, c1, t):
    a, b = to_rgb(c0), to_rgb(c1)
    return tuple(a[i] + (b[i] - a[i]) * t for i in range(3))


RAMP_LO = S.FEATURE_COLORS["theta"]        # P1 (#9AA7B3)
RAMP_HI = S.FEATURE_COLORS["ESM+theta"]    # P5 (#2E3742)
P_COLOR = {f"P{i+1}": _mix(RAMP_LO, RAMP_HI, i / 4.0) for i in range(5)}
# text ramp compressed to mid->dark slate (same rule as panel d) so the 5.8 pt
# P1..P5 labels stay legible while keeping the lightness ordering
P_TEXT = {f"P{i+1}": _mix(S.FEATURE_COLORS["ESM"], RAMP_HI, i / 4.0) for i in range(5)}
OBS = S.GREY                     # observed-score markers / intervals
LINK_MISS = S.GREY               # connectors, near-floor board (dashed)
LINK_HIT = S.INK                 # connectors, high-resolution board (solid)
PANEL_BG = "#F6F6F6"

TRUE_ORDER = ["P5", "P4", "P3", "P2", "P1"]        # top -> bottom
SCRAMBLED = ["P2", "P4", "P1", "P5", "P3"]         # observed order, near-floor benchmark
ROWS_Y = [28.0, 24.6, 21.2, 17.8, 14.4]


def rbox(ax, x0, y0, x1, y1, ec, fc="white", lw=0.5, r=0.8, z=1):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                boxstyle=f"round,pad=0,rounding_size={r}",
                                linewidth=lw, edgecolor=ec, facecolor=fc, zorder=z))


def mini_board(ax, x_true, x_obs_lab, x_bar0, bar_half, dots_x, obs_order,
               link_color, link_style):
    """One mini leaderboard: true-quality column, observed column, connecting lines."""
    for i, p in enumerate(TRUE_ORDER):
        y = ROWS_Y[i]
        ax.text(x_true, y - 0.75, p, fontsize=5.8, color=P_TEXT[p], ha="left",
                va="baseline", fontweight="bold")
        ax.scatter([x_true + 3.2], [y], s=6, color=P_COLOR[p], edgecolor="none",
                   zorder=4)

    for j, p in enumerate(obs_order):
        y = ROWS_Y[j]
        ax.text(x_obs_lab, y - 0.75, p, fontsize=5.8, color=P_TEXT[p], ha="left",
                va="baseline", fontweight="bold")
        xc = dots_x[j]
        ax.plot([xc - bar_half, xc + bar_half], [y, y], "-", color=OBS, lw=0.6,
                solid_capstyle="butt", zorder=3)
        ax.plot([xc - bar_half, xc - bar_half], [y - 0.6, y + 0.6], "-", color=OBS,
                lw=0.6, zorder=3)
        ax.plot([xc + bar_half, xc + bar_half], [y - 0.6, y + 0.6], "-", color=OBS,
                lw=0.6, zorder=3)
        ax.scatter([xc], [y], s=7, color=OBS, edgecolor="white", linewidths=0.3,
                   zorder=4)
        # connector from this predictor's true-quality row to its observed row
        i = TRUE_ORDER.index(p)
        ax.plot([x_true + 4.2, x_obs_lab - 0.8], [ROWS_Y[i], y], ls=link_style,
                color=link_color, lw=0.4, zorder=2)


def main() -> None:
    S.apply_rcparams()
    plt.rcParams["savefig.bbox"] = None

    fig, ax = S.panel(W_MM, H_MM)
    ax.set_position([0, 0, 1, 1])
    ax.set_xlim(0, W_MM)
    ax.set_ylim(0, H_MM)
    ax.axis("off")

    # divider between the two mini panels
    ax.plot([42.5, 42.5], [10.5, 32.0], "-", color=S.LIGHT_GREY, lw=0.5, zorder=1)

    # headers
    ax.text(1.0, 33.6, "Ceiling near chance (low resolution)", fontsize=5.8,
            fontweight="bold", color=S.INK, ha="left", va="baseline")
    ax.text(45.0, 33.6, "Ceiling above the floor (high resolution)", fontsize=5.8,
            fontweight="bold", color=S.INK, ha="left", va="baseline")

    # column sub-headers
    for x in (1.0, 45.0):
        ax.text(x, 30.6, "true quality", fontsize=5.5, color=S.GREY, ha="left",
                va="baseline")
    ax.text(19.0, 30.6, "observed score", fontsize=5.5, color=S.GREY, ha="left",
            va="baseline")
    ax.text(63.0, 30.6, "observed score", fontsize=5.5, color=S.GREY, ha="left",
            va="baseline")

    # left: overlapping intervals, order scrambled
    mini_board(ax, x_true=1.6, x_obs_lab=19.0, x_bar0=None, bar_half=3.6,
               dots_x=[30.2, 29.7, 29.2, 28.7, 28.2], obs_order=SCRAMBLED,
               link_color=LINK_MISS, link_style=(0, (2.2, 1.6)))
    # right: separated intervals, order recovered
    mini_board(ax, x_true=45.6, x_obs_lab=63.0, x_bar0=None, bar_half=1.1,
               dots_x=[78.4, 76.2, 74.0, 71.8, 69.6], obs_order=TRUE_ORDER,
               link_color=LINK_HIT, link_style="-")

    # caption strip
    rbox(ax, 1.0, 0.8, 84.0, 7.6, PANEL_BG, fc=PANEL_BG, lw=0.0, r=0.8, z=1)
    for y, t in ((5.0, "The truly better model need not rank higher when the measurement ceiling"),
                 (2.4, "is near chance; a higher-resolution benchmark recovers the true order.")):
        ax.text(42.5, y, t, fontsize=5.6, color=S.INK, ha="center", va="baseline")

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig5e_leaderboard"))


if __name__ == "__main__":
    main()
