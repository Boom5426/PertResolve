"""Figure 2e - DE fidelity declines as the required biological resolution increases.

Data-direct. Source: results/results_v4_exttheta.csv (same grid as 2c/2d/2f/2g) via
remote_data.exttheta(), in-house predictors. Reproducible gradient (row-mean over variants):
  direction agreement (sign concordance)       0.65
  DE-LFC Spearman (effect-size rank correlation) 0.25
  DE overlap (variant-specific gene identity)   0.15
These are DIFFERENT quantities (two fractions and one rank correlation), all bounded by
1 = perfect agreement with the measured response, which is what the x axis now states.
Message: predictors recover the program-level direction but progressively lose the finer
allele-specific signal.

(These replace the earlier manuscript values 78/46/28, which reproduced from no committed or
remote file; the .tex was updated to match.)

Nature Methods pass:
  * The former x label "Score (different quantities; see row labels)" pushed the
    reader to the row labels; the axis now names what is measured and its 1 = perfect
    anchor, and each row label still carries its own unit.
  * The former blue ramp (#8FB8DE / #5185C0 / #2C5A8F) reused TP53's gene hue for a
    non-gene variable. Colour encodes nothing here, so the bars are a single neutral
    grey and the coarse-to-fine ladder is carried by row order and the arrow.

Run:  python fig2e_de_gradient.py  ->  fig2e_de_gradient.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

HERE = os.path.dirname(os.path.abspath(__file__))

W_MM, H_MM = 90.5, 44.0
AX = (26.0, 10.5, 59.0, 29.0)   # left, bottom, width, height in mm

# (two-line y label, column, display) coarse -> fine.
# ``display`` is the statistic TYPE and drives both the bar style and the value
# format: proportions are filled bars labelled as percentages, the one rank
# correlation is an open bar labelled with rho, so the three quantities cannot be
# read as three points on a single proportion scale.
STEPS = [
    ("Direction agreement\n(sign, fraction)", "direction_agreement", "pct"),
    ("DE-LFC Spearman\n(rank correlation)", "DE_LFC_spearman", "corr"),
    ("DE overlap\n(gene set, fraction)", "DE_overlap", "pct"),
]


def canvas():
    fig = plt.figure(figsize=(W_MM * S.MM, H_MM * S.MM))
    bg = fig.add_axes([0, 0, 1, 1])
    bg.set_axis_off()
    bg.set_xlim(0, 1)
    bg.set_ylim(0, 1)
    bg.add_patch(Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="none",
                           zorder=-10))
    l, b, w, h = AX
    return fig, fig.add_axes([l / W_MM, b / H_MM, w / W_MM, h / H_MM])


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
    g = D.exttheta()
    inh = g[g.method.map(D.is_inhouse)]
    vals = [inh[c].mean() for _, c, _ in STEPS]

    fig, ax = canvas()
    ys = list(range(len(STEPS)))
    for y, (lab, col, disp), v in zip(ys, STEPS, vals):
        if disp == "pct":
            ax.barh(y, v, height=0.42, color=S.GREY, edgecolor="none", zorder=3)
            txt = f"{v*100:.0f}%"
        else:
            ax.barh(y, v, height=0.42, facecolor="white", edgecolor=S.GREY,
                    linewidth=0.6, zorder=3)
            txt = rf"$\rho$ = {v:.2f}"
        ax.text(v + 0.018, y, txt, va="center", ha="left", fontsize=6,
                color=S.INK)

    ax.set_yticks(ys)
    ax.set_yticklabels([lab for lab, *_ in STEPS], fontsize=5.6)
    ax.set_ylim(-0.55, len(STEPS) - 0.45)
    ax.invert_yaxis()  # coarse at top
    ax.set_xlim(0, 0.86)
    ax.set_xticks([0, 0.25, 0.5, 0.75])
    ax.set_xlabel("Agreement with measured response\n"
                  r"(fraction, filled; Spearman $\rho$, open; 1 = perfect)",
                  labelpad=1.5, linespacing=1.3)
    S.despine(ax, keep=("bottom",))
    ax.tick_params(length=2.2, axis="x", pad=1.5)
    ax.tick_params(length=0, axis="y", pad=1.5)

    # coarse -> fine resolution arrow on the right
    ax.annotate("", xy=(0.83, len(STEPS) - 0.7), xytext=(0.83, 0.7),
                arrowprops=dict(arrowstyle="->", lw=0.7, color=S.GREY))
    ax.text(0.855, len(STEPS) / 2 - 0.5, "increasing resolution", rotation=90,
            va="center", ha="left", fontsize=5.6, color=S.GREY)

    S.save(fig, os.path.join(HERE, "fig2e_de_gradient"))
    print("DE gradient (in-house exttheta):", [round(v, 3) for v in vals])


if __name__ == "__main__":
    main()
