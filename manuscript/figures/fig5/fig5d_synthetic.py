"""Figure 5d - Construction of synthetic predictors of known, graded quality (schematic).

Source file: none. Conceptual panel; no measured quantity is plotted. Profile shapes
are fixed illustrative literals. The interpolation shown is the one defined in the
manuscript Methods ("Benchmark resolution"):
    delta_hat(v, alpha) = (1 - alpha) * gene-mean + alpha * variant-specific delta,
so true predictor quality increases monotonically with alpha and the true ordering is
fixed by construction. Recovery probabilities are panels f and g, not here.
Numbers reproduced: none (the alpha tags 0/0.25/0.5/0.75/1 are the construction grid,
not results).
Drawn at 85 x 36 mm, i.e. its exact placement size in fig5_assemble.tex (scale 1.0).

Run:  python fig5d_synthetic.py  ->  fig5d_synthetic.pdf (+ .png)
Vector check:  pdfimages -list fig5d_synthetic.pdf | tail -n +3 | wc -l   # -> 0
"""
from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

W_MM, H_MM = 85.0, 36.0
# Sequential house slate ramp (nm_style.FEATURE_COLORS), light -> dark with true
# predictor quality. Gene hues are reserved for TP53/KRAS/GATA1/JAK1 in panels
# b, c, f, g, so this schematic borrows none of them; the ramp is monotone in
# lightness and therefore also survives greyscale. Panel e re-uses the SAME ramp
# for P1..P5, so d and e read as one two-step story.
GREY_P = S.FEATURE_COLORS["theta"]        # gene-mean end of the ramp (#9AA7B3)
BLUE = S.FEATURE_COLORS["ESM+theta"]      # variant-specific end of the ramp (#2E3742)
# Text uses a compressed version of the same ramp (mid slate -> dark slate) so that
# 5.5-6 pt labels keep >=5:1 contrast on white while still reading as one ramp.
TXT_LO = S.FEATURE_COLORS["ESM"]          # #5F6B76
TXT_HI = S.FEATURE_COLORS["ESM+theta"]    # #2E3742
PANEL_BG = "#F6F6F6"
CHIP_BG = "#EEF1F4"

ALPHAS = [0.0, 0.25, 0.5, 0.75, 1.0]
ALPHA_TAGS = ["0", "0.25", "0.5", "0.75", "1"]
PRED_TAGS = ["P1", "P2", "P3", "P4", "P5"]

# illustrative shapes (fixed literals, not data)
MEAN_SHAPE = [0.9, 0.5, 1.0, 0.4, 0.8, 0.6, 0.9, 0.5, 0.7, 0.6]
VAR_SHAPE = [1.6, -1.2, 0.4, 1.9, -0.6, 1.3, -1.7, 0.8, 1.5, -0.9]


def mix(c0, c1, t):
    a, b = to_rgb(c0), to_rgb(c1)
    return tuple(a[i] + (b[i] - a[i]) * t for i in range(3))


def rbox(ax, x0, y0, x1, y1, ec, fc="white", lw=0.5, r=0.9, z=2):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                boxstyle=f"round,pad=0,rounding_size={r}",
                                linewidth=lw, edgecolor=ec, facecolor=fc, zorder=z))


def profile(ax, x0, x1, y0, vals, color, amp, lw=0.8, ref=None):
    n = len(vals)
    step = (x1 - x0) / n
    scale = max(abs(v) for v in (ref if ref is not None else vals))
    ax.plot([x0, x1], [y0, y0], "-", color=S.LIGHT_GREY, lw=0.4, zorder=2)
    for i, v in enumerate(vals):
        x = x0 + step * (i + 0.5)
        ax.plot([x, x], [y0, y0 + amp * v / scale], "-", color=color, lw=lw,
                solid_capstyle="butt", zorder=3)


def main() -> None:
    S.apply_rcparams()
    plt.rcParams["savefig.bbox"] = None
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

    fig, ax = S.panel(W_MM, H_MM)
    ax.set_position([0, 0, 1, 1])
    ax.set_xlim(0, W_MM)
    ax.set_ylim(0, H_MM)
    ax.axis("off")

    # ---- left: the two reference profiles ---------------------------------
    ax.text(1.0, 31.4, r"gene-mean $\bar{\delta}_g$  ($\alpha$ = 0)", fontsize=5.8,
            color=TXT_LO, ha="left", va="baseline")
    profile(ax, 1.0, 23.0, 28.0, MEAN_SHAPE, GREY_P, amp=1.6, ref=VAR_SHAPE)

    ax.text(1.0, 21.8, r"variant-specific $\delta_v$  ($\alpha$ = 1)", fontsize=5.8,
            color=BLUE, ha="left", va="baseline")
    profile(ax, 1.0, 23.0, 18.4, VAR_SHAPE, BLUE, amp=2.4, ref=VAR_SHAPE)

    # ---- middle: the interpolation and the alpha series -------------------
    rbox(ax, 28.5, 27.6, 70.5, 33.4, S.LIGHT_GREY, fc=CHIP_BG, r=0.8)
    ax.text(49.5, 29.6,
            r"$\hat{\delta}(v,\alpha)=(1-\alpha)\,\bar{\delta}_g+\alpha\,\delta_v$",
            fontsize=7.0, color=S.INK, ha="center", va="baseline")

    # 5.5 pt grey sub-header: identical style to panel e's "true quality" /
    # "observed score" sub-headers, so d and e read as one matched pair.
    ax.text(28.5, 24.6, "synthetic predictor family", fontsize=5.5, color=S.GREY,
            ha="left", va="baseline")

    x0, wd, gap = 28.5, 6.6, 2.0
    for k, a in enumerate(ALPHAS):
        xs = x0 + k * (wd + gap)
        col = mix(GREY_P, BLUE, a)
        vals = [(1 - a) * m + a * v for m, v in zip(MEAN_SHAPE, VAR_SHAPE)]
        profile(ax, xs, xs + wd, 17.6, vals, col, amp=2.2, lw=0.7, ref=VAR_SHAPE)
        ax.text(xs + wd / 2, 13.2, ALPHA_TAGS[k], fontsize=5.8, color=S.INK,
                ha="center", va="baseline")
        ax.text(xs + wd / 2, 10.6, PRED_TAGS[k], fontsize=5.8,
                color=mix(TXT_LO, TXT_HI, a), ha="center", va="baseline",
                fontweight="bold")
    ax.text(27.4, 13.2, r"$\alpha$", fontsize=5.8, color=S.INK, ha="right",
            va="baseline")

    # ---- right: true-quality ladder ---------------------------------------
    ax.text(73.0, 31.4, "true quality", fontsize=5.5, color=S.GREY, ha="left",
            va="baseline")
    ax.add_patch(FancyArrowPatch((75.0, 12.4), (75.0, 29.4), arrowstyle="-|>",
                                 mutation_scale=4.0, linewidth=0.6, color=S.GREY,
                                 shrinkA=0, shrinkB=0, zorder=3))
    for k, a in enumerate(ALPHAS):
        y = 13.6 + k * 3.7
        col = mix(GREY_P, BLUE, a)
        ax.scatter([75.0], [y], s=9, color=col, edgecolor="white", linewidths=0.4,
                   zorder=5)
        ax.text(77.0, y - 0.8, PRED_TAGS[k], fontsize=5.8,
                color=mix(TXT_LO, TXT_HI, a), ha="left", va="baseline",
                fontweight="bold")

    # ---- caption strip ----------------------------------------------------
    rbox(ax, 1.0, 0.8, 84.0, 7.6, PANEL_BG, fc=PANEL_BG, lw=0.0, r=0.8, z=1)
    for y, t in ((5.0, "Interpolating between the gene-mean profile and each variant's measured"),
                 (2.4, "effect builds predictors whose true quality order is fixed by construction.")):
        ax.text(42.5, y, t, fontsize=5.6, color=S.INK, ha="center", va="baseline")

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig5d_synthetic"))


if __name__ == "__main__":
    main()
