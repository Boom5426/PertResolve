"""Fig. 2b — allele-discrimination forest plot.

The 18 feature-model heads are grouped by model family rather than ranked by
their observed PDS. Gene-mean and WT-null form a separate reference block.
Every 95% bootstrap interval overlaps the analytic chance level of 0.50.

This panel owns the shared method-label column and feature-space key used by
Fig. 2c. Run: python fig2b_pds_forest.py

Static QA contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42; outputs
.svg, .pdf, .png and .tiff at dpi=600. Final assembled figure target:
width_mm = 183.
"""
from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

HERE = os.path.dirname(os.path.abspath(__file__))

W_MM, H_MM = 58.0, 63.0
AX_BOTTOM_MM, AX_HEIGHT_MM = 8.0, 47.5
AX_LEFT_MM, AX_WIDTH_MM = 18.0, 38.5

LW = 0.85
MS = 2.8
CAP = 0.21

FEATURE_SEQUENCE = ("theta", "esm", "esm+theta")
FEAT_TEX = {"theta": r"$\theta$", "esm": "ESM",
            "esm+theta": r"ESM+$\theta$"}
METHOD_BAND = "#F5F7F8"
REFERENCE_BAND = "#FAF6EC"


def short_label(method: str) -> str:
    if method in D.REFS:
        return method
    head, _, feat = method.partition("-")
    return f"{head}  {FEAT_TEX[feat]}"


def row_layout() -> tuple[list[str], list[float]]:
    """Return bottom-to-top methods and aligned y positions for panels b/c."""
    top_to_bottom = [
        f"{head}-{feat}" for head in D.HEADS for feat in FEATURE_SEQUENCE
    ] + ["Gene-mean", "WT-null"]
    methods = list(reversed(top_to_bottom))
    # Uniform spacing maximizes final-size legibility; the reference block is
    # separated by its pale background rather than by compressing the rows.
    positions = [0.0, 1.0]
    positions += [2.0 + i + 0.12 * (i // 3) for i in range(18)]
    return methods, positions


ROW_YLIM = (-0.48, 20.08)


def row_order() -> list[str]:
    """Compatibility helper used by related figure scripts."""
    return row_layout()[0]


def color_of(method: str) -> str:
    if method in D.REFS:
        return S.FEATURE_COLORS["reference"]
    return S.FEATURE_COLORS[D.feature_space(method)]


def interval(ax, lo, hi, val, y, c):
    ax.plot([lo, hi], [y, y], color=c, lw=LW, solid_capstyle="butt", zorder=3)
    for xx in (lo, hi):
        ax.plot([xx, xx], [y - CAP, y + CAP], color=c, lw=LW,
                solid_capstyle="butt", zorder=3)
    ax.plot([val], [y], marker="o", ms=MS, mfc=c, mec="white", mew=0.4,
            ls="none", zorder=4)


def add_row_bands(ax) -> None:
    methods, ys = row_layout()
    pos = dict(zip(methods, ys))
    ax.axhspan(-0.48, 1.48, color=REFERENCE_BAND, zorder=0)
    for i, head in enumerate(D.HEADS):
        if i % 2:
            members = [f"{head}-{feat}" for feat in FEATURE_SEQUENCE]
            vals = [pos[m] for m in members]
            ax.axhspan(min(vals) - 0.48, max(vals) + 0.48,
                       color=METHOD_BAND, zorder=0)


def canvas():
    fig = plt.figure(figsize=(W_MM * S.MM, H_MM * S.MM))
    bg = fig.add_axes([0, 0, 1, 1])
    bg.set_axis_off()
    bg.set_xlim(0, 1)
    bg.set_ylim(0, 1)
    bg.add_patch(Rectangle((0, 0), 1, 1, facecolor="white",
                           edgecolor="none", zorder=-10))
    ax = fig.add_axes([AX_LEFT_MM / W_MM, AX_BOTTOM_MM / H_MM,
                       AX_WIDTH_MM / W_MM, AX_HEIGHT_MM / H_MM])
    return fig, ax


def main() -> None:
    S.apply_rcparams()
    df = D.definitive().set_index("method")
    methods, ys = row_layout()

    assert len(methods) == 20
    assert set(D.EXTERNAL).isdisjoint(methods)
    assert df.loc[methods, "crosses"].all()

    fig, ax = canvas()
    add_row_bands(ax)
    ax.axvline(0.50, color=S.INK, lw=0.6, ls=(0, (3, 2)), zorder=2)

    for method, y in zip(methods, ys):
        row = df.loc[method]
        interval(ax, row.ci_lo, row.ci_hi, row.PDS, y, color_of(method))

    ax.set_yticks(ys)
    ax.set_yticklabels([short_label(m) for m in methods], fontsize=5.3)
    for tick, method in zip(ax.get_yticklabels(), methods):
        tick.set_color(S.INK)
        if method in D.REFS:
            tick.set_style("italic")

    ax.set_ylim(*ROW_YLIM)
    ax.set_xlim(0.455, 0.550)
    ax.set_xticks([0.46, 0.50, 0.54])
    ax.set_xlabel(r"PDS$_{cos}$ (allele discrimination)", labelpad=1.5)
    ax.tick_params(axis="y", length=0, pad=1.3)
    ax.tick_params(axis="x", pad=1.5)
    S.despine(ax, keep=("left", "bottom"))
    ax.spines["left"].set_bounds(ROW_YLIM[0], max(ys) + 0.35)

    def dot(color, label):
        return Line2D([0], [0], marker="o", ms=MS, mfc=color, mec="white",
                      mew=0.4, ls="none", label=label)

    handles = [
        dot(S.FEATURE_COLORS["theta"], r"$\theta$"),
        dot(S.FEATURE_COLORS["ESM"], "ESM"),
        dot(S.FEATURE_COLORS["ESM+theta"], r"ESM+$\theta$"),
        dot(S.FEATURE_COLORS["reference"], "reference"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.018, 0.995),
               ncol=4, columnspacing=0.65, handlelength=0.45,
               handletextpad=0.25, borderpad=0, borderaxespad=0,
               fontsize=5.4, frameon=False)
    fig.text(0.625, 0.889, "chance", ha="center", va="center",
             fontsize=5.1, color=S.INK)
    fig.text(0.975, 0.925, "19/19 non-null CIs include 0.50",
             ha="right", va="center", fontsize=5.3, color=S.INK,
             fontstyle="italic")

    S.save(fig, os.path.join(HERE, "fig2b_pds_forest"), exact=True,
           formats=("pdf", "png", "svg", "tiff"))


if __name__ == "__main__":
    main()
