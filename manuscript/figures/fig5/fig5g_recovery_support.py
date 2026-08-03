"""Fig. 5g — controlled recovery by ceiling, depth and test-set support."""
from __future__ import annotations

W_MM, H_MM = 182.5, 43.0

from matplotlib.lines import Line2D

import fig5_common as S


def _draw_metric(ax, data, metric, marker, filled):
    for gene in S.GENE_ORDER:
        group = data[data["gene"] == gene].sort_values("depth_m")
        color = S.GENE_COLORS[gene]
        ax.plot(group["ceiling_pds"], group[metric], color=color, lw=0.8,
                alpha=0.48, zorder=1)
        for _, row in group.iterrows():
            face = color if filled else "white"
            ax.scatter(
                [row["ceiling_pds"]],
                [row[metric]],
                s=S.point_size(row["n_var"], minimum=10.0, maximum=40.0),
                marker=marker,
                facecolor=face,
                edgecolor=color,
                linewidths=0.85 if not filled else 0.45,
                zorder=4,
            )


def main() -> None:
    S.apply_style()
    data = S.read_csv("controlled_recovery.csv")

    fig = S.figure(W_MM, H_MM)
    grid = fig.add_gridspec(
        2, 1, left=0.13, right=0.79, bottom=0.20, top=0.95,
        hspace=0.20, height_ratios=[1, 1]
    )
    top = fig.add_subplot(grid[0, 0])
    bottom = fig.add_subplot(grid[1, 0], sharex=top)

    _draw_metric(top, data, "P_correct_order", "o", True)
    _draw_metric(bottom, data, "P_winner", "D", False)

    for ax in (top, bottom):
        ax.axvline(0.5, color=S.MID_GREY, lw=0.65, ls=(0, (3, 2)), zorder=0)
        ax.set_xlim(0.46, 0.98)
        ax.set_ylim(-0.04, 1.06)
        ax.set_yticks([0.0, 0.5, 1.0])
        S.despine(ax)
    top.tick_params(axis="x", labelbottom=False)
    top.set_ylabel("P(recover\nfull order)")
    bottom.set_ylabel("P(select\nbest)")
    bottom.set_xlabel("Replicate ceiling (PDS; gene × depth condition)")
    bottom.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9])

    # Direct labels occupy the otherwise unused right margin.
    fig.text(0.825, 0.735, "genes", fontsize=5.6, color=S.GREY,
             ha="left", va="top")
    for index, gene in enumerate(S.GENE_ORDER):
        y = 0.690 - index * 0.055
        fig.add_artist(Line2D([0.825, 0.855], [y, y], transform=fig.transFigure,
                              color=S.GENE_COLORS[gene], lw=1.4))
        fig.text(0.868, y, gene, fontsize=5.7, color=S.INK,
                 ha="left", va="center")

    fig.text(0.825, 0.435, "variants scored", fontsize=5.6, color=S.GREY,
             ha="left", va="top")
    for index, value in enumerate([10, 100, 250]):
        y = 0.383 - index * 0.070
        fig.add_artist(Line2D([0.84], [y], transform=fig.transFigure,
                              marker="o", linestyle="none",
                              markersize=(S.point_size(value, minimum=10.0,
                                                       maximum=40.0) ** 0.5) * 0.86,
                              markerfacecolor=S.MID_GREY,
                              markeredgecolor="white", markeredgewidth=0.4))
        fig.text(0.868, y, str(value), fontsize=5.5, color=S.INK,
                 ha="left", va="center")

    fig.text(0.825, 0.205, "depth m range", fontsize=5.6, color=S.GREY,
             ha="left", va="top")
    ranges = {"TP53": "25–50", "KRAS": "25–150",
              "GATA1": "25–250", "JAK1": "25–150"}
    for index, gene in enumerate(S.GENE_ORDER):
        y = 0.145 - index * 0.034
        fig.text(0.825, y, gene, fontsize=5.0, color=S.GREY,
                 ha="left", va="center")
        fig.text(0.930, y, ranges[gene], fontsize=5.0, color=S.INK,
                 ha="right", va="center")

    S.save(fig, "fig5g_recovery_support")


if __name__ == "__main__":
    main()
