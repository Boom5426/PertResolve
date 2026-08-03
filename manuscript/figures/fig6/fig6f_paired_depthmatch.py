"""Fig. 6f — paired native versus matched-50 detection-unrankable fractions."""
# Imported export contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42;
# outputs .svg, .pdf and .tiff at dpi=600 via fig6_common.save.
# Final assembled figure target: width_mm = 183.
from __future__ import annotations

W_MM, H_MM = 80.5, 37.0

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import fig6_common as S

ORDER = ["Replogle", "Adamson", "Norman"]
Y = {"Replogle": 2.0, "Adamson": 1.0, "Norman": 0.0}


def main() -> None:
    S.apply_style()
    data = pd.read_csv(S.DERIVED / "depthmatch_summary.csv").set_index("dataset")

    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0.29, 0.225, 0.67, 0.68])

    for dataset in ORDER:
        row = data.loc[dataset]
        y = Y[dataset]
        colour = S.DATASET_COLORS[dataset]
        native_y, matched_y = y + 0.09, y - 0.09
        native = float(row["native_pct"])
        matched = float(row["matched50_pct"])
        ax.plot(
            [native, matched],
            [native_y, matched_y],
            color=colour,
            lw=0.75,
            alpha=0.85,
            zorder=1,
        )
        ax.errorbar(
            native,
            native_y,
            xerr=np.array([[native - row["native_lo"]], [row["native_hi"] - native]]),
            fmt="o",
            ms=4.2,
            mfc=colour,
            mec="white",
            mew=0.45,
            ecolor=colour,
            elinewidth=0.75,
            capsize=1.5,
            capthick=0.65,
            zorder=4,
        )
        ax.errorbar(
            matched,
            matched_y,
            xerr=np.array([[matched - row["matched50_lo"]], [row["matched50_hi"] - matched]]),
            fmt="o",
            ms=4.5,
            mfc="white",
            mec=colour,
            mew=0.9,
            ecolor=colour,
            elinewidth=0.75,
            capsize=1.5,
            capthick=0.65,
            zorder=3,
        )
        annotation_x = max(float(row["native_hi"]), float(row["matched50_hi"])) + 2.2
        ax.text(
            annotation_x,
            y,
            f"Δ +{row['delta_pp']:.1f} pp",
            fontsize=5.6,
            color=S.INK,
            fontweight="bold",
            va="center",
            ha="left",
        )

    ax.set_xlim(0, 70)
    ax.set_ylim(-0.45, 2.45)
    ax.set_xticks([0, 20, 40, 60])
    ax.set_yticks([2, 1, 0])
    ax.set_yticklabels(
        [f"{dataset}  n={int(data.loc[dataset, 'paired_n']):,}" for dataset in ORDER]
    )
    for tick, dataset in zip(ax.get_yticklabels(), ORDER):
        tick.set_color(S.DATASET_COLORS[dataset])
        tick.set_fontweight("bold")
    ax.set_xlabel("Detection-unrankable perturbations (%)")
    S.despine(ax, keep=("bottom",))

    handles = [
        Line2D([0], [0], marker="o", linestyle="none", ms=4.0, mfc=S.DARK_SLATE, mec="white", mew=0.45, label="native depth"),
        Line2D([0], [0], marker="o", linestyle="none", ms=4.2, mfc="white", mec=S.DARK_SLATE, mew=0.8, label="matched 50 cells per half"),
    ]
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.29, 0.99),
        ncol=1,
        labelspacing=0.34,
        handletextpad=0.35,
        fontsize=5.5,
    )

    S.save(fig, "fig6f_paired_depthmatch")


if __name__ == "__main__":
    main()
