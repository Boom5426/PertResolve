"""Fig. 4a — allele-ranking robustness across generalization splits.

The previous 18 x 5 heatmap encoded small deviations around chance with nearly
white cells, so the central result required reading the colour bar. This
revision shows the 18 matched feature-model heads directly. Split-level means
and IQRs make the verdict visible without treating between-head dispersion as
an inferential confidence interval.

Static QA contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42; outputs
.svg, .pdf, .png and .tiff at dpi=600. Final assembled figure target:
width_mm = 183.
"""
from __future__ import annotations

import numpy as np
import matplotlib as mpl
from matplotlib.patches import Rectangle

import fig4_common as S


W_MM, H_MM = 89.0, 44.0
SPLITS = S.SPLIT_ORDER
POINT_COLOUR = "#8B949B"
LINE_COLOUR = "#D6DADD"
LOW_DEPTH_FACE = "#F7ECE8"


def summaries(matrix: np.ndarray):
    return {
        "mean": np.mean(matrix, axis=0),
        "median": np.median(matrix, axis=0),
        "q1": np.percentile(matrix, 25, axis=0),
        "q3": np.percentile(matrix, 75, axis=0),
        "minimum": np.min(matrix, axis=0),
        "maximum": np.max(matrix, axis=0),
    }


def main() -> None:
    S.apply_style()
    # Font keys deliberately not restated: fig4_common already resolves the
    # family through nm_style, and a local "Arial, Helvetica, DejaVu Sans" stack
    # silently renders in DejaVu Sans wherever Arial is absent.
    mpl.rcParams.update({
        "font.size": 5.0,
    })

    frame = S.breakdown()
    inhouse = frame[frame["method"].map(S.is_inhouse)]
    pivot = inhouse.groupby(["method", "split"])["PDS_cos"].mean().unstack()
    rows = S.ordered_methods()
    matrix = pivot.loc[rows, SPLITS].to_numpy(float)
    if matrix.shape != (18, 5) or not np.isfinite(matrix).all():
        raise ValueError(f"expected a finite 18 x 5 PDS matrix, observed {matrix.shape}")

    stats = summaries(matrix)
    four_near_chance = np.delete(stats["mean"], 3)
    max_non_low_depth_deviation = float(np.max(np.abs(four_near_chance - 0.50)))
    assert max_non_low_depth_deviation < 0.022
    assert stats["mean"][3] < 0.25

    fig = S.figure(W_MM, H_MM)
    bg = fig.add_axes([0, 0, 1, 1])
    bg.set_axis_off()
    bg.add_patch(Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="none"))
    ax = fig.add_axes([0.125, 0.255, 0.82, 0.625])

    x = np.arange(len(SPLITS), dtype=float)
    ax.axvspan(2.58, 3.42, facecolor=LOW_DEPTH_FACE, edgecolor="none", zorder=0)
    ax.axhline(0.50, color=S.INK, lw=0.65, ls=(0, (3, 2)), zorder=1)

    # Same feature-model head appears at every x position. Lines are kept very
    # light because the split-level distribution, not any single trajectory,
    # is the evidence of interest.
    for values in matrix:
        ax.plot(x, values, color=LINE_COLOUR, lw=0.45, alpha=0.72, zorder=2)

    offsets = np.linspace(-0.105, 0.105, matrix.shape[0])
    for column in range(matrix.shape[1]):
        point_colour = S.NEGATIVE if column == 3 else POINT_COLOUR
        ax.scatter(
            np.full(matrix.shape[0], x[column]) + offsets,
            matrix[:, column],
            s=9.5,
            facecolor=point_colour,
            edgecolor="white",
            linewidth=0.28,
            alpha=0.76,
            zorder=4,
        )
        mean_colour = S.NEGATIVE if column == 3 else S.PDS_BLUE
        ax.vlines(x[column], stats["q1"][column], stats["q3"][column],
                  color=mean_colour, lw=2.2, zorder=5)
        ax.hlines([stats["q1"][column], stats["q3"][column]],
                  x[column] - 0.075, x[column] + 0.075,
                  color=mean_colour, lw=1.0, zorder=5)
        ax.scatter([x[column]], [stats["mean"][column]], s=23, marker="D",
                   facecolor=mean_colour, edgecolor="white", linewidth=0.45,
                   zorder=6)

        label_y = 0.548 if column != 3 else stats["mean"][column] + 0.062
        ax.text(
            x[column], label_y, f"{stats['mean'][column]:.3f}",
            ha="center", va="center", fontsize=5.2, color=mean_colour,
            fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="none", pad=0.45),
            zorder=7,
        )

    ax.text(4.35, 0.507, "chance = 0.50", ha="right", va="bottom",
            fontsize=5.0, color=S.GREY,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.35), zorder=7)

    labels = [
        "Random\n4 genes",
        "Positional\n4 genes",
        "Mechanistic\n4 genes",
        "Low-depth\n2 genes",
        "Compatibility\n3 genes",
    ]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=5.2, linespacing=1.05)
    ax.get_xticklabels()[3].set_color(S.NEGATIVE)
    ax.get_xticklabels()[3].set_fontweight("bold")
    ax.set_xlim(-0.45, 4.45)
    ax.set_ylim(0.20, 0.565)
    ax.set_yticks([0.20, 0.30, 0.40, 0.50])
    ax.set_ylabel(r"PDS$_{cos}$  (allele discrimination)", labelpad=1.6)
    ax.tick_params(axis="x", length=0, pad=2.0)
    ax.tick_params(axis="y", pad=1.4)
    S.despine(ax, keep=("left", "bottom"))

    fig.text(0.535, 0.910,
             "points: 18 matched feature–model heads · diamond: mean · thick line: IQR",
             ha="center", va="top", fontsize=5.0, color=S.GREY)
    fig.text(
        0.535, 0.085,
        "Four split means lie within 0.021 of chance; low-depth falls to 0.244.",
        ha="center", va="center", fontsize=5.25, fontweight="bold", color=S.INK,
    )
    fig.text(
        0.535, 0.035,
        "Head means pool eligible held-out variants; gene composition differs by split.",
        ha="center", va="center", fontsize=5.0, color=S.GREY,
    )

    S.save(fig, "fig4a_pds_splits")


if __name__ == "__main__":
    main()
