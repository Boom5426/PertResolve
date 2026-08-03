"""Vector heatmap helpers shared by refined Fig. 4a and Fig. 4b."""
from __future__ import annotations

import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory

import fig4_common as S


def draw_method_heatmap(ax, matrix: np.ndarray) -> None:
    n_rows, n_cols = matrix.shape
    for row in range(n_rows):
        for column in range(n_cols):
            ax.add_patch(
                Rectangle(
                    (column - 0.5, row - 0.5),
                    1.0,
                    1.0,
                    facecolor=ax._fig3_cmap(ax._fig3_norm(matrix[row, column])),
                    edgecolor="white",
                    linewidth=0.45,
                    zorder=2,
                )
            )

    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(
        [
            f"{S.SPLIT_LABELS[split]}\n({S.SPLIT_GENE_COUNTS[split]} genes)"
            for split in S.SPLIT_ORDER
        ],
        rotation=55,
        ha="right",
        rotation_mode="anchor",
        fontsize=5.0,
        linespacing=1.05,
    )
    rows = S.ordered_methods()
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([S.head_name(method) for method in rows], fontsize=5.1)
    ax.tick_params(axis="both", length=0, pad=1.3)

    transform = blended_transform_factory(ax.transAxes, ax.transData)
    block_size = len(S.HEADS)
    for index, feature in enumerate(S.FEATURE_ORDER):
        ax.text(
            -0.29,
            index * block_size + (block_size - 1) / 2.0,
            S.FEATURE_LABELS[feature],
            transform=transform,
            rotation=90,
            ha="center",
            va="center",
            fontsize=6.2,
            color=S.INK,
        )
    for index in range(1, len(S.FEATURE_ORDER)):
        ax.axhline(index * block_size - 0.5, color=S.INK, lw=0.65, zorder=4)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)


def draw_vector_colorbar(cax, cmap, norm, vmin: float, vmax: float, ticks, labels) -> None:
    edges = np.linspace(vmin, vmax, 121)
    for low, high in zip(edges[:-1], edges[1:]):
        value = 0.5 * (low + high)
        cax.add_patch(
            Rectangle(
                (0.0, low),
                1.0,
                high - low,
                facecolor=cmap(norm(value)),
                edgecolor="none",
            )
        )
    cax.set_xlim(0.0, 1.0)
    cax.set_ylim(vmin, vmax)
    cax.set_xticks([])
    cax.set_yticks(ticks)
    cax.set_yticklabels(labels, fontsize=5.1)
    cax.yaxis.tick_right()
    cax.tick_params(axis="y", length=1.8, width=0.5, pad=1.3)
    for side in ("top", "right", "left", "bottom"):
        cax.spines[side].set_visible(True)
        cax.spines[side].set_linewidth(0.45)
        cax.spines[side].set_edgecolor(S.INK)
