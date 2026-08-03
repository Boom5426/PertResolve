"""Fig. 6b — retrospective same-source LODO ROC curves."""
# Imported export contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42;
# outputs .svg, .pdf and .tiff at dpi=600 via fig6_common.save.
# Final assembled figure target: width_mm = 183.
from __future__ import annotations

W_MM, H_MM = 65.5, 41.0

import pandas as pd

import fig6_common as S

DASHES = {
    "Replogle": (0, ()),
    "Norman": (0, (4.0, 1.5)),
    "Adamson": (0, (1.1, 1.3)),
    "GATA1": (0, (4.5, 1.3, 1.0, 1.3)),
    "JAK1": (0, (2.2, 1.1, 0.8, 1.1, 0.8, 1.1)),
}
ORDER = ["Replogle", "Norman", "Adamson", "GATA1", "JAK1"]


def main() -> None:
    S.apply_style()
    curves = pd.read_csv(S.DERIVED / "retrospective_roc.csv")
    predictions = pd.read_csv(S.DERIVED / "retrospective_predictions.csv")

    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0.11, 0.255, 0.49, 0.66])
    key = fig.add_axes([0.655, 0.24, 0.33, 0.66])
    key.set_xlim(0, 1)
    key.set_ylim(0, 1)
    key.axis("off")

    ax.plot([0, 1], [0, 1], ls=(0, (2.0, 2.0)), lw=0.6, color=S.MID_GREY, zorder=1)
    for dataset in ORDER:
        group = curves[curves["dataset"] == dataset]
        ax.plot(
            group["fpr"],
            group["tpr"],
            color=S.DATASET_COLORS[dataset],
            ls=DASHES[dataset],
            lw=1.15,
            solid_capstyle="round",
            dash_capstyle="butt",
            zorder=3,
        )

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_xlabel("False-positive rate")
    ax.set_ylabel("True-positive rate")
    ax.set_aspect("equal", adjustable="box")
    S.despine(ax)

    # Two header lines, not two headers on one line: at 29 mm the column cannot
    # carry "Held-out dataset" and "rankable / unrankable" side by side, and
    # stacking them on one line put the second over the first data row.
    key.text(0.0, 0.99, "Held-out dataset", fontsize=5.5, color=S.GREY,
             fontweight="bold", va="top")
    key.text(0.98, 0.90, "rankable / unrankable", fontsize=5.5, color=S.GREY,
             fontweight="bold", ha="right", va="top")
    for row, dataset in enumerate(ORDER):
        y = 0.68 - row * 0.155
        subset = predictions[predictions["dataset"] == dataset]
        n_pos = int(subset["n_pos"].iloc[0])
        n_neg = int(subset["n_neg"].iloc[0])
        key.plot(
            [0.00, 0.22],
            [y, y],
            color=S.DATASET_COLORS[dataset],
            ls=DASHES[dataset],
            lw=1.2,
            solid_capstyle="round",
            dash_capstyle="butt",
        )
        # 6.2 pt name plus a 5.8 pt count needs more than the 24 mm key column
        # gives; the name drops a size so the two clear each other.
        key.text(0.28, y, dataset, fontsize=5.8, va="center", color=S.DATASET_COLORS[dataset], fontweight="bold")
        key.text(0.98, y, f"{n_pos:,} / {n_neg:,}", fontsize=5.8, va="center", ha="right", color=S.INK)

    S.save(fig, "fig6b_retrospective_roc")


if __name__ == "__main__":
    main()
