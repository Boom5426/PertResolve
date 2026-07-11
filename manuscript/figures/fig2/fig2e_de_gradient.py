"""Figure 2e - DE fidelity declines as the required biological resolution increases.

Data-direct. Source: results/results_v4_exttheta.csv (same grid as 2c/2d/2f/2g) via
remote_data.exttheta(), in-house predictors. Reproducible gradient (row-mean over variants):
  direction agreement (sign concordance)       0.65
  DE-LFC Spearman (effect-size ranking)         0.26
  DE overlap (variant-specific gene identity)   0.15
These are DIFFERENT quantities (two fractions and one rank correlation), labelled as such
and read as a coarse->fine ladder, not a single comparable axis. Message: predictors recover
the program-level direction but progressively lose the finer allele-specific signal.

(These replace the earlier manuscript values 78/46/28, which reproduced from no committed or
remote file; the .tex was updated to match.)

Run:  python fig2e_de_gradient.py  ->  fig2e_de_gradient.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

# (two-line y label, column, display) coarse -> fine
STEPS = [
    ("Direction agreement\n(sign, fraction)", "direction_agreement", "pct"),
    ("DE-LFC Spearman\n(effect-size ranking)", "DE_LFC_spearman", "corr"),
    ("DE overlap\n(variant-specific genes)", "DE_overlap", "pct"),
]
SHADES = ["#8FB8DE", "#5185C0", "#2C5A8F"]  # coarse (light) -> fine (dark)


def main() -> None:
    S.apply_rcparams()
    g = D.exttheta()
    inh = g[g.method.map(D.is_inhouse)]
    vals = [inh[c].mean() for _, c, _ in STEPS]

    fig, ax = S.panel(66, 44)
    ys = list(range(len(STEPS)))
    for y, (lab, col, disp), v, sh in zip(ys, STEPS, vals, SHADES):
        ax.barh(y, v, height=0.6, color=sh, edgecolor="none", zorder=3)
        txt = f"{v*100:.0f}%" if disp == "pct" else f"{v:.2f}"
        ax.text(v + 0.012, y, txt, va="center", ha="left", fontsize=7,
                color=S.INK, fontweight="bold")

    ax.set_yticks(ys)
    ax.set_yticklabels([lab for lab, *_ in STEPS], fontsize=5.6)
    ax.set_ylim(-0.6, len(STEPS) - 0.4)
    ax.invert_yaxis()  # coarse at top
    ax.set_xlim(0, 0.82)
    ax.set_xticks([0, 0.25, 0.5, 0.75])
    ax.set_xlabel("Score (different quantities; see row labels)")
    S.despine(ax, keep=("bottom",))
    ax.tick_params(length=2.2, axis="x")
    ax.tick_params(length=0, axis="y")

    # coarse -> fine resolution arrow on the right
    ax.annotate("", xy=(0.79, len(STEPS) - 0.75), xytext=(0.79, 0.75),
                arrowprops=dict(arrowstyle="->", lw=0.8, color=S.GREY))
    ax.text(0.815, len(STEPS) / 2 - 0.5, "increasing resolution", rotation=90,
            va="center", ha="left", fontsize=5, color=S.GREY)

    S.save(fig, "fig2e_de_gradient")
    print("DE gradient (in-house exttheta):", [round(v, 3) for v in vals])


if __name__ == "__main__":
    main()
