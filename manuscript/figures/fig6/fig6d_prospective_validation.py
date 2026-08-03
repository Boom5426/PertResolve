"""Fig. 6d — prospective pilot-to-future validation and calibration."""
# Imported export contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42;
# outputs .svg, .pdf and .tiff at dpi=600 via fig6_common.save.
# Final assembled figure target: width_mm = 183.
from __future__ import annotations

W_MM, H_MM = 80.5, 43.0

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import fig6_common as S

DATASETS = ["Replogle", "Norman", "Adamson"]
Y = {"Replogle": 2.0, "Norman": 1.0, "Adamson": 0.0}
SERIES = [
    ("Learned effect", 25, -0.12, "o", "white", 3.7),
    ("Learned effect", 50, 0.12, "o", "filled", 4.2),
    ("Training-free SNR", 50, 0.00, "D", "white", 3.4),
]


def main() -> None:
    S.apply_style()
    auc = pd.read_csv(S.DERIVED / "prospective_auc.csv")
    calibration = pd.read_csv(S.DERIVED / "prospective_calibration.csv")

    fig = S.figure(W_MM, H_MM)
    forest = fig.add_axes([0.155, 0.235, 0.42, 0.61])
    calib = fig.add_axes([0.685, 0.245, 0.225, 0.57])

    forest.axvline(0.5, ls=(0, (2.0, 2.0)), lw=0.6, color=S.MID_GREY, zorder=1)
    for dataset in DATASETS:
        colour = S.DATASET_COLORS[dataset]
        for method, cells, offset, marker, face_mode, size in SERIES:
            row = auc[
                (auc["dataset"] == dataset)
                & (auc["method"] == method)
                & (auc["pilot_cells"] == cells)
            ].iloc[0]
            value, low, high = float(row["auroc"]), float(row["lo"]), float(row["hi"])
            face = colour if face_mode == "filled" else "white"
            forest.errorbar(
                value,
                Y[dataset] + offset,
                xerr=np.array([[value - low], [high - value]]),
                fmt=marker,
                ms=size,
                mfc=face,
                mec=colour,
                mew=0.8,
                ecolor=colour,
                elinewidth=0.75,
                capsize=1.5,
                capthick=0.65,
                zorder=4 if face_mode == "filled" else 3,
            )

    forest.set_xlim(0.48, 1.02)
    forest.set_ylim(-0.45, 2.45)
    forest.set_xticks([0.5, 0.75, 1.0])
    forest.set_yticks([2, 1, 0])
    forest.set_yticklabels(DATASETS)
    for tick, dataset in zip(forest.get_yticklabels(), DATASETS):
        tick.set_color(S.DATASET_COLORS[dataset])
        tick.set_fontweight("bold")
    forest.set_xlabel("Prospective AUROC")
    S.despine(forest, keep=("bottom",))

    handles = [
        Line2D([0], [0], marker="o", linestyle="none", ms=3.8, mfc="white", mec=S.DARK_SLATE, mew=0.8, label="25-cell learned effect"),
        Line2D([0], [0], marker="o", linestyle="none", ms=4.1, mfc=S.DARK_SLATE, mec=S.DARK_SLATE, label="50-cell learned effect"),
        Line2D([0], [0], marker="D", linestyle="none", ms=3.3, mfc="white", mec=S.DARK_SLATE, mew=0.8, label="50-cell training-free SNR"),
    ]
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.105, 0.985),
        ncol=1,
        labelspacing=0.34,
        handletextpad=0.35,
        fontsize=5.4,
    )

    calib.plot([0, 1], [0, 1], ls=(0, (2.0, 2.0)), lw=0.65, color=S.MID_GREY, zorder=1)
    for _, row in calibration.iterrows():
        calib.errorbar(
            row["mean_predicted"],
            row["observed_rate"],
            yerr=np.array(
                [[row["observed_rate"] - row["lo"]], [row["hi"] - row["observed_rate"]]]
            ),
            fmt="o",
            ms=4.0,
            mfc=S.DARK_SLATE,
            mec="white",
            mew=0.45,
            ecolor=S.DARK_SLATE,
            elinewidth=0.75,
            capsize=1.5,
            capthick=0.65,
            zorder=3,
        )
    calib.set_xlim(-0.03, 1.03)
    calib.set_ylim(-0.03, 1.03)
    calib.set_xticks([0, 0.5, 1.0])
    calib.set_yticks([0, 0.5, 1.0])
    calib.set_xlabel("Mean predicted probability")
    calib.set_ylabel("Observed rankable fraction")
    calib.set_aspect("equal", adjustable="box")
    S.despine(calib)
    calib.text(0.0, 1.20, "Calibration", transform=calib.transAxes, fontsize=6.2, fontweight="bold", va="bottom")
    # Right-aligned on the axis: left-aligned this ran 12.9 mm past the 89 mm
    # canvas once the panel narrowed.
    calib.text(1.0, 1.045, "50-cell learned model", transform=calib.transAxes,
               fontsize=5.1, color=S.GREY, va="bottom", ha="right")


    S.save(fig, "fig6d_prospective_validation")


if __name__ == "__main__":
    main()
