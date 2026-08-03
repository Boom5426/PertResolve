"""Fig. 2e — biological fidelity across progressively finer resolutions.

The panel is a biological-resolution ladder, not a common-scale bar chart.
Direction agreement, DE-LFC Spearman correlation and DE-set overlap retain
their native definitions and nulls. Values and observed ranges are calculated
across all 18 in-house feature-model heads.

Static QA contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42; outputs
.svg, .pdf, .png and .tiff at dpi=600. Final assembled figure target:
width_mm = 183.
"""
from __future__ import annotations

import os
import sys

import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

HERE = os.path.dirname(os.path.abspath(__file__))
W_MM, H_MM = 89.5, 44.0

TEAL = S.RESOLUTION
TEAL_FACE = S.RESOLUTION_LIGHT
DARK = S.FEATURE_COLORS["ESM+theta"]
MID = S.FEATURE_COLORS["ESM"]
PALE = "#F4F6F7"
CENTERS = [15.1, 45.25, 75.4]

STAGES = [
    {
        "title": "Program direction",
        "subtitle": "sign agreement",
        "column": "direction_agreement",
        "value": lambda v: f"{v * 100:.0f}%",
        "range": lambda a, b: f"{a * 100:.0f}–{b * 100:.0f}%",
        "null": "chance = 50%",
        "status": "retained",
    },
    {
        "title": "Effect-size order",
        "subtitle": "DE-LFC Spearman",
        "column": "DE_LFC_spearman",
        "value": lambda v: rf"$\rho$ = {v:.2f}",
        "range": lambda a, b: rf"$\rho$ {a:.2f}–{b:.2f}",
        "null": r"null $\rho$ = 0",
        "status": "weakly ordered",
    },
    {
        "title": "Responding-gene identity",
        "subtitle": "DE overlap",
        "column": "DE_overlap",
        "value": lambda v: f"{v * 100:.0f}%",
        "range": lambda a, b: f"{a * 100:.0f}–{b * 100:.0f}%",
        "null": "top-50 overlap\nset-size dependent",
        "status": "largely lost",
    },
]


def draw_sign_glyph(ax, cx, cy):
    """Small observed/predicted sign patterns; conceptual, not extra data."""
    measured = np.array([1, -1, 1, 1, -1, 1, -1])
    predicted = np.array([1, -1, 1, -1, -1, 1, 1])
    xs = np.linspace(cx - 5.2, cx + 6.2, len(measured))
    ax.text(cx - 8.7, cy + 1.45, "obs.", ha="left", va="center",
            fontsize=5.0, color=S.GREY)
    ax.text(cx - 8.7, cy - 1.45, "pred.", ha="left", va="center",
            fontsize=5.0, color=S.GREY)
    for x, obs, pred in zip(xs, measured, predicted):
        match = obs == pred
        colour = TEAL if match else S.LIGHT_GREY
        ax.scatter([x], [cy + 1.45], marker="^" if obs > 0 else "v", s=9,
                   facecolor=colour, edgecolor="white", linewidth=0.25, zorder=4)
        ax.scatter([x], [cy - 1.45], marker="^" if pred > 0 else "v", s=9,
                   facecolor=colour, edgecolor="white", linewidth=0.25, zorder=4)
        ax.plot([x, x], [cy - 0.9, cy + 0.9], color=colour, lw=0.45, zorder=2)


def draw_rank_glyph(ax, cx, cy):
    """Two partially discordant effect-size rankings."""
    left_order = [0, 1, 2, 3, 4]
    right_order = [1, 0, 3, 2, 4]
    ys = np.linspace(cy + 3.3, cy - 3.3, 5)
    x_left, x_right = cx - 8.0, cx + 3.0
    ax.text(x_left + 2.3, cy + 5.2, "observed", ha="center", va="center",
            fontsize=5.0, color=S.GREY)
    ax.text(x_right + 2.3, cy + 5.2, "predicted", ha="center", va="center",
            fontsize=5.0, color=S.GREY)
    widths = [4.6, 4.0, 3.4, 2.8, 2.2]
    for item in range(5):
        yl = ys[left_order.index(item)]
        yr = ys[right_order.index(item)]
        ax.plot([x_left + 4.8, x_right - 0.3], [yl, yr], color=S.LIGHT_GREY,
                lw=0.45, zorder=1)
        colour = TEAL if item in (0, 4) else MID
        ax.add_patch(Rectangle((x_left, yl - 0.38), widths[item], 0.76,
                               facecolor=colour, edgecolor="none", alpha=0.78))
        ax.add_patch(Rectangle((x_right, yr - 0.38), widths[item], 0.76,
                               facecolor=colour, edgecolor="none", alpha=0.78))


def draw_overlap_glyph(ax, cx, cy):
    """Predicted and observed DE sets with a small shared region."""
    ax.add_patch(Circle((cx - 2.2, cy), 4.3, facecolor=PALE,
                        edgecolor=MID, linewidth=0.55, zorder=1))
    ax.add_patch(Circle((cx + 2.2, cy), 4.3, facecolor=PALE,
                        edgecolor=MID, linewidth=0.55, zorder=1))
    ax.add_patch(Ellipse((cx, cy), 1.7, 6.3, facecolor=TEAL_FACE,
                         edgecolor=TEAL, linewidth=0.55, zorder=2))
    ax.text(cx - 6.1, cy + 5.15, "predicted\nDE set", ha="center", va="center",
            fontsize=5.0, color=S.GREY, linespacing=1.0)
    ax.text(cx + 6.1, cy + 5.15, "observed\nDE set", ha="center", va="center",
            fontsize=5.0, color=S.GREY, linespacing=1.0)


def main() -> None:
    S.apply_rcparams()
    frame = D.exttheta()
    inhouse = frame[frame.method.map(D.is_inhouse)]
    per_method = inhouse.groupby("method")[[s["column"] for s in STAGES]].mean()
    assert per_method.shape[0] == 18
    if not np.isfinite(per_method.to_numpy(float)).all():
        raise ValueError("DE-fidelity summaries contain non-finite values")

    fig, ax = S.panel(W_MM, H_MM)
    ax.set_position([0, 0, 1, 1])
    ax.set_xlim(0, W_MM)
    ax.set_ylim(0, H_MM)
    ax.set_axis_off()
    ax.add_patch(Rectangle((0, 0), W_MM, H_MM, facecolor="white",
                           edgecolor="none", zorder=-10))

    ax.text(1.6, 41.1, "coarser gene-shared signal", ha="left", va="center",
            fontsize=5.2, color=S.GREY)
    ax.text(W_MM - 1.6, 41.1, "finer allele-specific resolution", ha="right",
            va="center", fontsize=5.2, color=S.GREY)
    ax.add_patch(FancyArrowPatch((24.0, 41.1), (63.8, 41.1),
                                arrowstyle="-|>", mutation_scale=4.5,
                                lw=0.65, color=S.GREY))

    for x in (30.15, 60.30):
        ax.plot([x, x], [7.2, 38.8], color=S.LIGHT_GREY, lw=0.5)

    for cx, stage in zip(CENTERS, STAGES):
        values = per_method[stage["column"]].to_numpy(float)
        mean = float(values.mean())
        lo, hi = float(values.min()), float(values.max())
        ax.text(cx, 37.7, stage["title"], ha="center", va="center",
                fontsize=5.8, fontweight="bold", color=S.INK)
        ax.text(cx, 35.3, stage["subtitle"], ha="center", va="center",
                fontsize=5.0, color=S.GREY)
        ax.text(cx, 31.3, stage["value"](mean), ha="center", va="center",
                fontsize=9.0, fontweight="bold", color=S.INK)

    draw_sign_glyph(ax, CENTERS[0], 23.2)
    draw_rank_glyph(ax, CENTERS[1], 22.7)
    draw_overlap_glyph(ax, CENTERS[2], 22.7)

    for cx, stage in zip(CENTERS, STAGES):
        values = per_method[stage["column"]].to_numpy(float)
        lo, hi = float(values.min()), float(values.max())
        ax.text(cx, 15.0, f"18 heads · range {stage['range'](lo, hi)}",
                ha="center", va="center", fontsize=5.0, color=S.INK)
        ax.text(cx, 12.1, stage["null"], ha="center", va="center",
                fontsize=5.0, color=S.GREY, linespacing=1.0)
        ax.text(cx, 9.0, stage["status"], ha="center", va="center",
                fontsize=5.1, color=TEAL if stage["status"] == "retained" else DARK,
                fontweight="bold")

    ax.text(W_MM / 2, 5.4,
            "Gene-shared direction is retained; finer allele-specific gene identity is largely lost.",
            ha="center", va="center", fontsize=5.25, color=S.INK,
            fontweight="bold")
    ax.text(W_MM / 2, 2.0,
            "Native metrics have different nulls; compare biological resolution, not raw magnitude.",
            ha="center", va="center", fontsize=5.0, color=S.GREY,
            fontstyle="italic")

    S.save(fig, os.path.join(HERE, "fig2e_de_gradient"), exact=True,
           formats=("pdf", "png", "svg", "tiff"))


if __name__ == "__main__":
    main()
