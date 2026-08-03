"""Figure 1e — model scoring after the measurement-identification gate.

The panel makes two levels explicit: single-cell identification establishes an
allele-benchmarkable ground truth in Fig. 1a, whereas PDS is a model score used
inside that validated benchmark. Direction recovery and DE-program fidelity
remain supporting axes; allele discrimination is the primary model endpoint.

Static QA contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42; outputs
.svg, .pdf, .png and .tiff at dpi=600. Final assembled figure target:
width_mm = 183.
"""
from __future__ import annotations

import os
import sys

import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

from pathlib import Path

HERE = Path(__file__).resolve().parent

W_MM, H_MM = 78.0, 37.0
ID_EDGE = "#688896"
ID_FACE = "#EAF1F4"
PALE_GREY = "#F5F6F7"


def rounded_box(ax, x, y, width, height, *, face="white", edge=S.LIGHT_GREY,
                lw=0.55, radius=0.8, zorder=0):
    patch = FancyBboxPatch(
        (x, y), width, height,
        boxstyle=f"round,pad=0.08,rounding_size={radius}",
        facecolor=face, edgecolor=edge, linewidth=lw, zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def direction_icon(ax, cx):
    """Observed and predicted profiles share direction, but not allele identity."""
    x = np.linspace(cx - 7.5, cx + 6.1, 70)
    phase = np.linspace(0, 2.3 * np.pi, 70)
    y_pred = 25.4 + 1.45 * np.sin(phase + 0.10)
    y_obs = 22.5 + 1.28 * np.sin(phase)
    ax.plot(x, y_pred, color=S.INK, lw=0.7)
    ax.plot(x, y_obs, color=S.GREY, lw=0.7, ls=(0, (2.0, 1.35)))
    label_box = {"facecolor": "white", "edgecolor": "none", "pad": 0.05,
                 "alpha": 0.94}
    ax.text(cx + 7.25, y_pred[-1], "pred.", fontsize=5.0, ha="center",
            va="center", bbox=label_box)
    ax.text(cx + 7.25, y_obs[-1], "obs.", fontsize=5.0, ha="center",
            va="center", color=S.GREY, bbox=label_box)


def pds_icon(ax, cx):
    """The prediction retrieves its own observed allele above its siblings.

    The whole group sits 1.9 mm further left than it was first drawn. The card
    is 32.2 mm wide and the row must hold the prediction box, the arrow, the
    longest bar and the label "1  own observed"; laid out from the old origin
    that label crossed the card's right border. Shifting the group rebalances
    the two insets to about 2 mm each instead of 4.1 and -0.1.
    """
    rounded_box(ax, cx - 13.9, 23.2, 4.7, 3.5, face="white", edge=S.INK,
                lw=0.5, radius=0.6, zorder=2)
    ax.text(cx - 11.55, 24.95, r"$\hat{\delta}_v$", fontsize=5.5,
            ha="center", va="center", zorder=3)
    ax.add_patch(
        FancyArrowPatch(
            (cx - 8.8, 24.95), (cx - 6.4, 24.95), arrowstyle="-|>",
            mutation_scale=4.5, lw=0.5, color=S.GREY,
        )
    )
    rows = [
        (26.7, 8.1, S.RESOLUTION, "1  own observed"),
        (24.4, 5.8, S.LIGHT_GREY, "2  sibling A"),
        (22.1, 4.5, S.LIGHT_GREY, "3  sibling B"),
    ]
    for y, width, colour, label in rows:
        ax.add_patch(Rectangle((cx - 6.6, y - 0.62), width, 1.24,
                               facecolor=colour, edgecolor="none", zorder=1))
        ax.text(cx + 2.7, y, label, fontsize=5.0,
                color=S.RESOLUTION if colour == S.RESOLUTION else S.GREY,
                ha="left", va="center")


def de_icon(ax, cx):
    """Predicted and observed changed-gene lists share a response program."""
    for x0, label, colour in [(cx - 8.2, "pred.", S.INK),
                              (cx + 1.1, "obs.", S.GREY)]:
        ax.text(x0 + 3.3, 26.5, label, fontsize=5.0, color=colour,
                ha="center", va="bottom")
        for i in range(4):
            y = 24.7 - i * 1.75
            shared = (x0 < cx and i in (0, 2)) or (x0 > cx and i in (1, 2))
            ax.add_patch(
                Rectangle(
                    (x0, y - 0.52), 6.6, 1.04,
                    facecolor=S.RESOLUTION if shared else S.LIGHT_GREY,
                    edgecolor="none",
                )
            )


def main() -> None:
    S.apply_rcparams()
    fig, ax = S.panel(W_MM, H_MM)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_xlim(0, W_MM)
    ax.set_ylim(0, H_MM)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(Rectangle((0, 0), W_MM, H_MM, facecolor="white",
                           edgecolor="none", zorder=-10))

    # Explicit hand-off from Fig. 1a: identification validates the target;
    # the three cards below score the model prediction, not the measurement.
    rounded_box(ax, 0.9, 34.0, 21.7, 2.15, face=ID_FACE, edge=ID_EDGE,
                lw=0.55, radius=0.5)
    ax.text(11.75, 35.08, "IDENTIFICATION PASS (a)", fontsize=5.0,
            color=ID_EDGE, fontweight="bold", ha="center", va="center")
    ax.add_patch(
        FancyArrowPatch(
            (23.1, 35.08), (26.2, 35.08), arrowstyle="-|>",
            mutation_scale=4.6, lw=0.55, color=S.GREY,
        )
    )
    ax.text(27.0, 35.08, "SCORE THE HELD-OUT PREDICTION", fontsize=5.15,
            color=S.INK, fontweight="bold", ha="left", va="center")

    card_specs = [
        (0.8, 21.3, "white", S.LIGHT_GREY),
        (22.9, 32.2, S.RESOLUTION_LIGHT, S.RESOLUTION),
        (56.2, 21.0, "white", S.LIGHT_GREY),
    ]
    for x, width, face, edge in card_specs:
        rounded_box(ax, x, 4.25, width, 28.75, face=face, edge=edge,
                    lw=0.55, radius=0.8)

    centers = [11.45, 39.0, 66.7]
    headers = ["Direction recovery", "Allele discrimination", "DE-program fidelity"]
    for i, (cx, header) in enumerate(zip(centers, headers)):
        ax.text(cx, 31.35, header, fontsize=5.8 if i != 1 else 6.2,
                fontweight="bold", color=S.INK, ha="center", va="center")
        ax.text(cx, 29.0, "PRIMARY MODEL SCORE" if i == 1 else "SUPPORTING AXIS",
                fontsize=5.0, fontweight="bold",
                color=S.RESOLUTION if i == 1 else S.GREY,
                ha="center", va="center")

    direction_icon(ax, centers[0])
    pds_icon(ax, centers[1])
    de_icon(ax, centers[2])

    ax.text(centers[0], 18.25, "Pearson-δ · δ-cosine", fontsize=5.25,
            ha="center")
    ax.text(centers[0], 13.1, "Is the shared response", fontsize=5.0,
            color=S.GREY, ha="center")
    ax.text(centers[0], 10.65, "direction recovered?", fontsize=5.0,
            color=S.GREY, ha="center")
    ax.text(centers[1], 18.5, "PDS", fontsize=6.2, color=S.INK,
            fontweight="bold", ha="center")
    ax.text(centers[1], 16.2, "tie-aware percentile rank",
            fontsize=5.0, ha="center")
    ax.text(centers[1], 14.55, "of its own observed allele",
            fontsize=5.0, ha="center")
    rounded_box(ax, 26.0, 10.9, 26.0, 2.75, face="white",
                edge=S.RESOLUTION, lw=0.45, radius=0.42)
    ax.text(centers[1], 12.28, "chance = 0.50 · higher is better", fontsize=5.0,
            color=S.INK, ha="center", va="center", fontweight="bold")
    ax.text(centers[1], 8.65, r"Does $\hat{\delta}_v$ rank its own observed allele",
            fontsize=5.0, color=S.GREY, ha="center")
    ax.text(centers[1], 6.55, "above sibling alleles?", fontsize=5.0,
            color=S.GREY, ha="center")

    ax.text(centers[2], 18.2, "DE overlap", fontsize=5.2, ha="center")
    ax.text(centers[2], 15.8, "DE-LFC rank", fontsize=5.0, ha="center")
    ax.text(centers[2], 13.4, "direction agreement", fontsize=5.0, ha="center")
    ax.text(centers[2], 9.4, "Do changed-gene", fontsize=5.0,
            color=S.GREY, ha="center")
    ax.text(centers[2], 7.3, "programs agree?", fontsize=5.0,
            color=S.GREY, ha="center")

    ax.text(
        W_MM / 2, 1.7,
        "Measurement identification validates the benchmark; PDS evaluates the model within it.",
        fontsize=5.0, color=S.INK, fontweight="bold", ha="center", va="center",
    )

    S.save(fig, HERE / "fig1e_eval_axes", exact=True,
           formats=("pdf", "png", "svg", "tiff"))


if __name__ == "__main__":
    main()
