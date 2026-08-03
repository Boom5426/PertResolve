"""Fig. 6e — exact native-depth effect-size landscape for allele datasets."""
# Imported export contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42;
# outputs .svg, .pdf and .tiff at dpi=600 via fig6_common.save.
# Final assembled figure target: width_mm = 183.
from __future__ import annotations

W_MM, H_MM = 80, 37.0

import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import fig6_common as S


def main() -> None:
    S.apply_style()
    data = pd.read_csv(S.DERIVED / "allele_native_rankability.csv")
    valid = (
        np.isfinite(data["effect_size"])
        & np.isfinite(data["n_cells"])
        & (data["effect_size"] > 0)
        & (data["n_cells"] > 0)
    )
    if not bool(valid.all()):
        raise ValueError("Panel e requires finite, strictly positive effect sizes and depths.")
    cap_count = int((data["n_cells"] == 300).sum())

    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0.175, 0.235, 0.785, 0.615])
    ax.axhspan(286, 315, color=S.PALE_GREY, zorder=0)
    ax.axhline(300, ls=(0, (1.4, 1.5)), lw=0.6, color=S.MID_GREY, zorder=1)

    for gene in S.GENE_ORDER:
        group = data[data["gene"] == gene]
        colour = S.DATASET_COLORS[gene]
        unrankable = group[~group["rankable"].astype(bool)]
        rankable = group[group["rankable"].astype(bool)]
        ax.scatter(
            unrankable["effect_size"],
            unrankable["n_cells"],
            s=6.0,
            facecolors="none",
            edgecolors=colour,
            linewidths=0.35,
            alpha=0.58,
            zorder=2,
        )
        ax.scatter(
            rankable["effect_size"],
            rankable["n_cells"],
            s=16.0,
            facecolors=colour,
            edgecolors="white",
            linewidths=0.45,
            alpha=0.98,
            zorder=4,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.62, 17.0)
    ax.set_ylim(55, 480)
    ax.set_xticks([1, 3, 10])
    ax.set_xticklabels(["1", "3", "10"])
    ax.set_yticks([100, 200, 300])
    ax.set_yticklabels(["100", "200", "300"])
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.yaxis.set_minor_locator(mticker.NullLocator())
    ax.set_xlabel(r"Effect size,  $\Vert\delta\Vert_2$")
    ax.set_ylabel("Cells scored per perturbation\n(deepest split-half rung)")
    S.despine(ax)
    ax.spines["left"].set_bounds(55, 300)
    ax.text(
        0.66,
        321,
        f"analysis cap · {cap_count}/{len(data)} perturbations",
        fontsize=5.2,
        color=S.GREY,
        ha="left",
        va="bottom",
    )

    # Compact, non-overlapping gene key in the deliberately empty header band.
    x_positions = [0.02, 0.27, 0.52, 0.77]
    for x, gene in zip(x_positions, S.GENE_ORDER):
        group = data[data["gene"] == gene]
        count = int(group["rankable"].sum())
        total = len(group)
        ax.text(
            x,
            1.06,
            f"{gene}  {count}/{total}",
            transform=ax.transAxes,
            color=S.DATASET_COLORS[gene],
            fontsize=5.7,
            fontweight="bold",
            ha="left",
            va="bottom",
            clip_on=False,
        )

    state_handles = [
        Line2D([0], [0], marker="o", linestyle="none", ms=3.4, mfc=S.GREY, mec="white", mew=0.4, label="rankable"),
        Line2D([0], [0], marker="o", linestyle="none", ms=3.0, mfc="white", mec=S.GREY, mew=0.5, label="detection-unrankable"),
    ]
    ax.legend(
        handles=state_handles,
        loc="lower left",
        bbox_to_anchor=(0.0, 0.0),
        handletextpad=0.3,
        labelspacing=0.25,
        borderpad=0.15,
        fontsize=5.4,
    )

    S.save(fig, "fig6e_effect_depth_landscape")


if __name__ == "__main__":
    main()
