"""Figure 5a - The oracle predictor defined (conceptual schematic, drawn at final size).

Source file: none. This panel is conceptual and plots no measured quantity; every
profile shape is a fixed illustrative literal, not a result. Per-gene ceiling values
live in panel b, so no number appears here by design.
Numbers reproduced: none (deliberate; see fig5a_prompt.md "Do NOT include").
Message: the oracle prediction for a held-out variant is an independent second
split-half measurement of that same variant, scored through the identical PDS
harness, so its score is the discrimination attainable at the available depth.
Drawn at 70 x 44 mm, i.e. its exact placement size in fig5_assemble.tex (scale 1.0),
so all type is natively 5.5-6 pt on the page.

Run:  python fig5a_oracle_def.py  ->  fig5a_oracle_def.pdf (+ .png)
Vector check:  pdfimages -list fig5a_oracle_def.pdf | tail -n +3 | wc -l   # -> 0
"""
from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

# Drawing canvas (coordinate units, mm-like) and the physical figure size. The
# canvas is drawn into a 64 x 40.2 mm figure, i.e. exactly its placement size in
# fig5_assemble.tex (scale 1.0), so every point size below is the on-page size.
W_MM, H_MM = 70.0, 44.0
FIG_W, FIG_H = 64.0, 41.1
# Neutral slate/grey pair, NOT gene hues: in this figure blue/orange/purple/green
# are reserved for TP53/KRAS/GATA1/JAK1 (panels b, c, f, g). Half 1 -> measured
# target (grey), half 2 -> oracle prediction (dark slate); the lightness step also
# survives greyscale.
MEAS = S.GREY                          # measured target / half 1
ORAC = S.FEATURE_COLORS["ESM+theta"]   # oracle prediction / half 2  (#2E3742)
PANEL_BG = "#F6F6F6"
ORAC_BG = "#EEF1F4"

# illustrative pseudobulk shapes (fixed literals, not data)
SHAPE_A = [1.7, -0.9, 2.2, 0.6, -1.8, 1.1, -0.4, 2.0, -1.3, 0.8, 1.5, -2.1, 0.5, 1.9]
SHAPE_B = [1.5, -1.1, 2.0, 0.9, -1.6, 1.3, -0.6, 1.8, -1.5, 0.6, 1.7, -1.9, 0.7, 1.6]
# two half-cluster dot layouts (fixed literals)
DOTS_1 = [(0.10, 0.72), (0.30, 0.28), (0.48, 0.80), (0.66, 0.36), (0.84, 0.66),
          (0.20, 0.16), (0.58, 0.10), (0.92, 0.20), (0.38, 0.55)]
DOTS_2 = [(0.14, 0.30), (0.34, 0.74), (0.52, 0.22), (0.70, 0.68), (0.88, 0.34),
          (0.24, 0.86), (0.62, 0.88), (0.94, 0.78), (0.44, 0.50)]


def rbox(ax, x0, y0, x1, y1, ec, fc="white", lw=0.5, r=0.9, z=2):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                boxstyle=f"round,pad=0,rounding_size={r}",
                                linewidth=lw, edgecolor=ec, facecolor=fc, zorder=z))


def arrow(ax, p0, p1, color, lw=0.6, rad=0.0, z=4):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=4.0,
                                 linewidth=lw, color=color, zorder=z,
                                 shrinkA=0, shrinkB=0,
                                 connectionstyle=f"arc3,rad={rad}"))


def profile(ax, x0, x1, y0, vals, color, amp=2.3, lw=0.9):
    """Signed mini bar profile around a zero line (illustrative)."""
    n = len(vals)
    step = (x1 - x0) / n
    m = max(abs(v) for v in vals)
    ax.plot([x0, x1], [y0, y0], "-", color=S.LIGHT_GREY, lw=0.4, zorder=2)
    for i, v in enumerate(vals):
        x = x0 + step * (i + 0.5)
        ax.plot([x, x], [y0, y0 + amp * v / m], "-", color=color, lw=lw,
                solid_capstyle="butt", zorder=3)


def dots(ax, x0, x1, y0, y1, layout, color):
    xs = [x0 + (x1 - x0) * u for u, _ in layout]
    ys = [y0 + (y1 - y0) * v for _, v in layout]
    ax.scatter(xs, ys, s=3.0, color=color, edgecolor="none", zorder=3)


def main() -> None:
    S.apply_rcparams()
    plt.rcParams["savefig.bbox"] = None   # exact canvas -> scale 1.0 in the composite
    # Keep mathtext in the same Arial-metric sans as the rest of the figure. Without
    # this, matplotlib renders $\delta$ / $\alpha$ from DejaVu Sans, so the panel would
    # ship two different typefaces. nm_style is shared and must not be edited, so the
    # override is local to this panel.
    plt.rcParams.update({
        "mathtext.fontset": "custom",
        "mathtext.rm": "Liberation Sans",
        "mathtext.it": "Liberation Sans:italic",
        "mathtext.bf": "Liberation Sans:bold",
    })

    fig, ax = S.panel(FIG_W, FIG_H)
    ax.set_position([0, 0, 1, 1])
    ax.set_xlim(0, W_MM)
    ax.set_ylim(0, H_MM)
    ax.axis("off")

    # ---- step headers -----------------------------------------------------
    # Header 3 is left-shifted and set one step smaller so it cannot overrun the
    # right edge of the 70-unit canvas at the 64 mm placement width.
    for x, t in ((1.0, "1  Split the cells"),
                 (21.0, "2  Two pseudobulks"),
                 (43.0, "3  Score in one harness")):
        ax.text(x, 41.0, t, fontsize=5.6, fontweight="bold", color=S.INK,
                ha="left", va="baseline")

    # ---- step 1: one variant's cells, split in two ------------------------
    rbox(ax, 1.5, 23.0, 18.5, 36.0, S.LIGHT_GREY)
    ax.text(10.0, 37.0, "variant v", fontsize=6.0, color=S.INK, ha="center",
            va="baseline", fontweight="bold")
    ax.plot([2.2, 17.8], [29.5, 29.5], ls=(0, (2, 1.6)), lw=0.5, color=S.GREY, zorder=3)
    ax.text(3.0, 32.4, "half 1", fontsize=5.5, color=MEAS, ha="left", va="center")
    ax.text(3.0, 26.2, "half 2", fontsize=5.5, color=ORAC, ha="left", va="center")
    dots(ax, 8.6, 17.4, 30.6, 35.0, DOTS_1, MEAS)
    dots(ax, 8.6, 17.4, 24.0, 28.4, DOTS_2, ORAC)

    arrow(ax, (18.9, 32.4), (20.7, 34.6), MEAS)
    arrow(ax, (18.9, 26.2), (20.7, 22.6), ORAC)

    # ---- step 2: two independent pseudobulks ------------------------------
    ax.text(21.0, 36.8, r"measured target $\delta(v)$", fontsize=5.8, color=MEAS,
            ha="left", va="baseline")
    profile(ax, 21.0, 38.0, 33.4, SHAPE_A, MEAS)

    ax.text(21.0, 27.4, "oracle prediction", fontsize=5.8, color=ORAC,
            ha="left", va="baseline")
    ax.text(21.0, 25.2, "(2nd measurement)", fontsize=5.8, color=ORAC,
            ha="left", va="baseline")
    profile(ax, 21.0, 38.0, 21.4, SHAPE_B, ORAC)

    arrow(ax, (38.8, 33.4), (44.6, 31.4), MEAS)
    arrow(ax, (38.8, 21.4), (44.6, 27.0), ORAC)

    # ---- step 3: identical harness -> ceiling -----------------------------
    rbox(ax, 45.0, 24.5, 69.0, 34.0, S.INK)
    ax.text(57.0, 31.0, "PDS harness", fontsize=6.0, fontweight="bold", color=S.INK,
            ha="center", va="baseline")
    ax.text(57.0, 28.3, "ranked against other", fontsize=5.5, color=S.GREY,
            ha="center", va="baseline")
    ax.text(57.0, 26.1, "held-out variants", fontsize=5.5, color=S.GREY,
            ha="center", va="baseline")

    arrow(ax, (57.0, 24.3), (57.0, 20.4), S.INK)

    rbox(ax, 45.0, 12.5, 69.0, 20.0, ORAC, fc=ORAC_BG, lw=0.8)
    ax.text(57.0, 16.9, "oracle discrimination", fontsize=6.0, color=ORAC,
            ha="center", va="baseline", fontweight="bold")
    ax.text(57.0, 14.3, "ceiling", fontsize=6.0, color=ORAC, ha="center",
            va="baseline", fontweight="bold")

    # ---- caption strip ----------------------------------------------------
    rbox(ax, 1.0, 1.2, 69.0, 9.6, PANEL_BG, fc=PANEL_BG, lw=0.0, r=0.8, z=1)
    for y, t in ((7.0, "The oracle prediction is a genuine second measurement"),
                 (4.8, "of the same variant, so its score is the discrimination"),
                 (2.6, "attainable at this depth.")):
        ax.text(35.0, y, t, fontsize=5.6, color=S.INK, ha="center", va="baseline")

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig5a_oracle_def"))


if __name__ == "__main__":
    main()
