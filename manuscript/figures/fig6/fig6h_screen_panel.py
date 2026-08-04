"""Fig. 6h — detection is common across published screens; identification is not.

Six perturbation screens run through the frozen criterion, three of which took no part in
developing it. The distance of each point below the diagonal is the gap between being able
to tell a perturbation from its control and being able to tell it from its closest
competitor, which is the gap the whole study is about.

Traces to results/canonical/resolution_panel.csv, produced under the pre-registration in
docs/PREREG_RESOLUTION_PANEL_v1.md.
"""
# Imported export contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42;
# outputs .svg, .pdf and .tiff at dpi=600 via fig6_common.save.
# Final assembled figure target: width_mm = 183.
from __future__ import annotations

W_MM, H_MM = 64.5, 38.0

import numpy as np

import fig6_common as S

# Where each label sits relative to its point, chosen so no two collide at this size.
LABELS = {
    "mcfarland":    ("McFarland",  -4.5, 4.0, "right"),
    "tahoe100m":    ("Tahoe-100M", -4.5, 3.0, "right"),
    "norman2019":   ("Norman",      5.0, 2.5, "left"),
    "replogle":     ("Replogle",    5.0, -3.0, "left"),
    "adamson2016":  ("Adamson",     3.0, 7.5, "left"),
    "vcc_training": ("VCC",        -3.0, 7.5, "right"),
}


def main() -> None:
    S.apply_style()
    panel = S.read_csv("resolution_panel.csv")

    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0.155, 0.235, 0.80, 0.62])

    # Below this line a screen detects more perturbations than it can tell apart. Every
    # screen measured here sits below it, most of them far below.
    ax.plot([0, 100], [0, 100], color=S.MID_GREY, lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax.text(52, 57, "equal", fontsize=5.0, color=S.GREY, rotation=45,
            ha="center", va="center")

    x = panel["detectable_fraction"] * 100.0
    y = panel["identifiable_fraction"] * 100.0
    sizes = np.clip(np.sqrt(panel["n_perturbations"]) * 2.6, 9.0, 30.0)
    ax.scatter(x, y, s=sizes, facecolor=S.BLUE_PALE, edgecolor=S.DARK_SLATE,
               linewidths=0.8, zorder=4)

    for _, row in panel.iterrows():
        label, dx, dy, ha = LABELS[row["dataset_name"]]
        ax.text(row["detectable_fraction"] * 100.0 + dx,
                row["identifiable_fraction"] * 100.0 + dy,
                label, fontsize=5.3, color=S.DARK_SLATE, ha=ha, va="center", zorder=5)

    ax.set_xlim(0, 100)
    ax.set_ylim(-3, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlabel("detectable (%)", labelpad=1.5)
    ax.set_ylabel("identifiable (%)", labelpad=2.0)
    S.despine(ax)

    # Anchored to the canvas edge, not the axes edge: at this width the axes inset alone left
    # too little room and the note ran 3.7 mm past the right side.
    fig.text(0.022, 0.985,
             "Six screens, frozen criterion at 50 cells per group; area, perturbations",
             fontsize=5.2, color=S.GREY, ha="left", va="top")
    S.save(fig, "fig6h_screen_panel")


if __name__ == "__main__":
    main()
