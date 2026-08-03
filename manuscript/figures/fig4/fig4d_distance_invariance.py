"""Fig. 4d — paired changes in PDS under alternative distance geometries."""
from __future__ import annotations

import numpy as np

import fig4_common as S

# Static QA contract: Arial/Helvetica sans-serif; svg.fonttype="none";
# pdf.fonttype=42; exports .svg, .pdf, .png and .tiff at dpi=600.
# Final assembled figure target: width_mm = 183.
W_MM, H_MM = 73.0, 37.0


def main() -> None:
    S.apply_style()
    df = S.breakdown()
    inhouse = df[df["method"].map(S.is_inhouse)]
    per_method = inhouse.groupby("method")[["PDS_cos", "PDS_L1", "PDS_L2"]].mean()
    differences = {
        "L1": (per_method["PDS_L1"] - per_method["PDS_cos"]).sort_values(),
        "L2": (per_method["PDS_L2"] - per_method["PDS_cos"]).sort_values(),
    }

    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0.25, 0.24, 0.70, 0.615])
    ax.axvline(0.0, color=S.INK, lw=0.7, ls=(0, (3, 2)), zorder=1)
    y_positions = {"L1": 1.0, "L2": 0.0}
    for label, values in differences.items():
        values = values.to_numpy(float)
        y = y_positions[label]
        offsets = S.deterministic_offsets(len(values), 0.12)
        ax.scatter(
            values,
            y + offsets,
            s=8,
            facecolor=S.MID_GREY,
            edgecolor="white",
            linewidth=0.3,
            zorder=3,
        )
        q1, median, q3 = np.quantile(values, [0.25, 0.50, 0.75])
        ax.plot([q1, q3], [y, y], color=S.SLATE, lw=2.1, zorder=4)
        ax.plot(
            median,
            y,
            marker="D",
            ms=4.1,
            mfc=S.PDS_BLUE,
            mec="white",
            mew=0.5,
            zorder=5,
        )
        raw_median = float(
            per_method["PDS_L1" if label == "L1" else "PDS_L2"].median()
        )
        ax.text(
            0.019,
            y,
            f"raw median {raw_median:.3f}",
            ha="right",
            va="center",
            fontsize=5.1,
            color=S.GREY,
        )

    maximum = max(np.abs(values.to_numpy(float)).max() for values in differences.values())
    ax.set_xlim(-0.020, 0.020)
    ax.set_ylim(-0.42, 1.42)
    ax.set_xticks([-0.02, -0.01, 0.00, 0.01, 0.02])
    ax.set_xticklabels(["−0.02", "−0.01", "0", "+0.01", "+0.02"])
    ax.set_yticks([1.0, 0.0])
    ax.set_yticklabels(["L1 − cosine", "L2 − cosine"])
    ax.set_xlabel("Within-head change in PDS")
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", pad=1.3)
    S.despine(ax, keep=("bottom",))
    ax.text(
        0.02,
        0.97,
        f"max |change| = {maximum:.3f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.2,
        color=S.INK,
    )
    fig.text(
        0.55,
        0.895,
        "points: 18 heads · diamond: median · thick line: IQR",
        ha="center",
        va="top",
        fontsize=5.0,
        color=S.GREY,
    )
    S.save(fig, "fig4d_distance_invariance")


if __name__ == "__main__":
    main()
