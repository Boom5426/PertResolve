"""Fig. 5b — native-depth replicate ceiling versus all in-house model heads."""
from __future__ import annotations

W_MM, H_MM = 78.0, 42.0

from matplotlib.lines import Line2D

import fig5_common as S


def main() -> None:
    S.apply_style()
    oracle = S.read_csv("oracle_ceiling.csv").query("scope != 'ALL'").set_index("scope")
    heads = S.per_gene_head_means()

    fig, ax = S.panel(W_MM, H_MM)
    fig.subplots_adjust(left=0.205, right=0.975, bottom=0.29, top=0.95)

    y_map = {gene: 3 - index for index, gene in enumerate(S.GENE_ORDER)}
    ax.axvline(0.5, color=S.MID_GREY, lw=0.65, ls=(0, (3, 2)), zorder=0)
    ax.text(0.5, 3.52, "chance", fontsize=5.4, color=S.GREY,
            ha="center", va="bottom")

    for gene in S.GENE_ORDER:
        row = oracle.loc[gene]
        values = heads.loc[heads["gene"] == gene, "PDS_cos"].to_numpy()
        y = y_map[gene]
        color = S.GENE_COLORS[gene]
        offsets = S.deterministic_offsets(len(values), 0.12)
        model_y = y - 0.10 + offsets
        best = float(values.max())
        ceiling = float(row["PDS_oracle"])

        ax.plot([best, ceiling], [y - 0.10, y + 0.13], color=S.LIGHT_GREY,
                lw=0.75, zorder=1)
        ax.scatter(values, model_y, s=9.5, color=color, alpha=0.42,
                   edgecolor="none", zorder=2)
        ax.scatter([best], [y - 0.10], s=30, marker="D", facecolor="white",
                   edgecolor=color, linewidths=0.9, zorder=5)
        ax.errorbar(
            [ceiling],
            [y + 0.13],
            xerr=[[ceiling - float(row["ci_lo"])], [float(row["ci_hi"]) - ceiling]],
            fmt="o",
            markersize=4.7,
            color=color,
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.45,
            ecolor=color,
            elinewidth=1.0,
            capsize=1.6,
            capthick=0.7,
            zorder=6,
        )

        if gene == "JAK1":
            gap = ceiling - best
            ax.text((ceiling + best) / 2.0, y + 0.29, f"gap = {gap:.3f}",
                    fontsize=5.5, color=S.GREY, ha="center", va="bottom")

    labels = [
        f"{gene}  n={int(oracle.loc[gene, 'n_var_scored'])}"
        for gene in S.GENE_ORDER
    ]
    ax.set_yticks([y_map[gene] for gene in S.GENE_ORDER], labels)
    ax.set_xlim(0.425, 0.865)
    ax.set_ylim(-0.55, 3.65)
    ax.set_xticks([0.5, 0.6, 0.7, 0.8])
    ax.set_xlabel("PDS (native-depth per-gene evaluation)")
    ax.tick_params(axis="y", length=0)
    S.despine(ax)

    handles = [
        Line2D([0], [0], marker="o", color=S.INK, markerfacecolor=S.INK,
               markeredgecolor="white", markersize=4.2, lw=0.8,
               label="replicate ceiling (95% CI)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=S.MID_GREY,
               markeredgecolor="none", alpha=0.55, markersize=3.3,
               label="18 model heads"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor="white",
               markeredgecolor=S.INK, markeredgewidth=0.8, markersize=4.0,
               label="maximum of 18"),
    ]
    # Two rows, not three across: at 78 mm the third entry ran 6 mm past the
    # panel. The row above the plot is the only free band, so it stacks.
    # Top-right corner, stacked: the axes occupy the lower two thirds of the
    # panel, so this band is free and the key never crosses the data.
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.995, 0.995),
               ncol=1, labelspacing=0.32, handletextpad=0.45)

    S.save(fig, "fig5b_ceiling_vs_models")


if __name__ == "__main__":
    main()

