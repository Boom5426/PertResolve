"""Fig. 4b — transcriptional direction recovery across five splits."""
from __future__ import annotations

import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize

import fig4_common as S
from fig4_heatmap import draw_method_heatmap, draw_vector_colorbar

# Static QA contract: Arial/Helvetica sans-serif; svg.fonttype="none";
# pdf.fonttype=42; exports .svg, .pdf, .png and .tiff at dpi=600.
# Final assembled figure target: width_mm = 183.
W_MM, H_MM = 89.0, 44.0


def main() -> None:
    S.apply_style()
    df = S.breakdown()
    inhouse = df[df["method"].map(S.is_inhouse)]
    pivot = inhouse.groupby(["method", "split"])["pearson_delta"].mean().unstack()
    rows = S.ordered_methods()
    matrix = pivot.loc[rows, S.SPLIT_ORDER].to_numpy(float)
    if np.isnan(matrix).any():
        raise ValueError("Fig. 4b matrix contains missing cells")

    cmap = LinearSegmentedColormap.from_list(
        "direction_score", ["#FFFFFF", "#DDE6EE", "#9FB6CA", S.DIRECTION_BLUE]
    )
    norm = Normalize(vmin=0.0, vmax=0.80)

    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0.23, 0.28, 0.57, 0.605])
    ax._fig3_cmap = cmap
    ax._fig3_norm = norm
    draw_method_heatmap(ax, matrix)

    cax = fig.add_axes([0.84, 0.30, 0.035, 0.50])
    ticks = [0.0, 0.2, 0.4, 0.6, 0.8]
    draw_vector_colorbar(
        cax,
        cmap,
        norm,
        0.0,
        0.80,
        ticks,
        ["0", "0.2", "0.4", "0.6", "0.8"],
    )
    fig.text(0.84, 0.82, r"Pearson-$\delta$", ha="left", va="bottom", fontsize=5.4)
    fig.text(
        0.52,
        0.915,
        "18 feature–model heads · independent scale with zero as baseline",
        ha="center",
        va="top",
        fontsize=5.1,
        color=S.GREY,
    )
    fig.text(
        0.23,
        0.035,
        "single-draw breakdown · cell: mean over held-out variants",
        ha="left",
        va="bottom",
        fontsize=5.0,
        color=S.GREY,
    )
    S.save(fig, "fig4b_direction_splits")


if __name__ == "__main__":
    main()
