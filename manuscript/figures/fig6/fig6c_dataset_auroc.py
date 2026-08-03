"""Fig. 6c — dataset-level retrospective LODO AUROC forest."""
# Imported export contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42;
# outputs .svg, .pdf and .tiff at dpi=600 via fig6_common.save.
# Final assembled figure target: width_mm = 183.
from __future__ import annotations

W_MM, H_MM = 80, 43.0

import numpy as np
import pandas as pd

import fig6_common as S

ROWS = [
    ("Replogle", 7.0),
    ("Norman", 6.0),
    ("Adamson", 5.0),
    ("GATA1", 3.2),
    ("JAK1", 2.2),
    ("TP53", 1.2),
    ("KRAS", 0.2),
]


def main() -> None:
    S.apply_style()
    data = pd.read_csv(S.RAW / "rankability_predictor_honest.csv")
    data = data[data["feature_set"] == "effect_size"].set_index("held_out")

    fig = S.figure(W_MM, H_MM)
    # Bottom raised from 0.14 to leave room for the JAK1 caveat below the
    # axis: placed at the old offset that note fell entirely outside the
    # 72 mm canvas (y -4.12 to -0.05 mm) and was clipped away completely.
    ax = fig.add_axes([0.40, 0.335, 0.56, 0.625])
    ax.axvline(0.5, ls=(0, (2.0, 2.0)), lw=0.6, color=S.MID_GREY, zorder=1)

    for dataset, y in ROWS:
        row = data.loc[dataset]
        n = int(row["n_perturbations"])
        n_pos = int(round(n * float(row["rankable_rate"])))
        n_neg = n - n_pos
        label = f"{dataset}   {n_pos:,}+ / {n_neg:,}−"
        ax.text(
            0.46,
            y,
            label,
            ha="right",
            va="center",
            fontsize=5.9,
            color=S.DATASET_COLORS[dataset],
            fontweight="bold" if bool(row["evaluable"]) else "normal",
            clip_on=False,
        )
        if bool(row["evaluable"]):
            value = float(row["auroc"])
            low = float(row["auroc_lo"])
            high = float(row["auroc_hi"])
            ax.errorbar(
                value,
                y,
                xerr=np.array([[value - low], [high - value]]),
                fmt="o",
                ms=4.0,
                mfc=S.DATASET_COLORS[dataset],
                mec="white",
                mew=0.45,
                ecolor=S.DATASET_COLORS[dataset],
                elinewidth=0.8,
                capsize=1.7,
                capthick=0.7,
                zorder=3,
            )
            ax.text(
                min(high + 0.018, 1.005),
                y,
                f"{value:.2f}",
                ha="left" if high + 0.018 <= 1.005 else "right",
                va="center",
                fontsize=5.4,
                color=S.INK,
            )
        else:
            ax.text(0.53, y, "NE · single class", ha="left", va="center", fontsize=5.5, color=S.GREY)

    ax.text(0.46, 7.72, "EXTERNAL ATLASES", ha="right", va="center", fontsize=5.1, color=S.GREY, fontweight="bold", clip_on=False)
    ax.text(0.46, 3.92, "ALLELE GENES", ha="right", va="center", fontsize=5.1, color=S.GREY, fontweight="bold", clip_on=False)
    ax.axhline(4.25, color=S.LIGHT_GREY, lw=0.55, zorder=0)
    ax.set_xlim(0.48, 1.03)
    ax.set_ylim(-0.35, 7.95)
    ax.set_xticks([0.5, 0.75, 1.0])
    ax.set_yticks([])
    ax.set_xlabel("Held-out AUROC")
    S.despine(ax, keep=("bottom",))
    fig.text(
        0.40,
        0.035,
        "JAK1: complete separation (18+ / 2−); the 1.00 estimate\nis based on two negatives.",
        fontsize=5.3,
        color=S.GREY,
        ha="left",
        va="bottom",
    )

    S.save(fig, "fig6c_dataset_auroc")


if __name__ == "__main__":
    main()
