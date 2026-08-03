"""Fig. 3d — pair-level and nearest-sibling allele identification."""
from __future__ import annotations

from matplotlib.lines import Line2D

import fig3_common as S

# Publication contract inherited from fig3_common: Arial/Helvetica sans-serif;
# svg.fonttype="none"; pdf.fonttype=42; exports .svg, .pdf, .png and .tiff at dpi=600.
# Final assembled figure target: width_mm = 183.
W_MM, H_MM = 100.0, 43.0


def main() -> None:
    S.apply_style()
    df = S.read_csv("pairwise_resolvability.csv").set_index("gene")
    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0.16, 0.24, 0.80, 0.60])

    y_positions = list(range(len(S.GENE_ORDER)))[::-1]
    dy = 0.20
    for y, gene in zip(y_positions, S.GENE_ORDER):
        color = S.GENE_COLORS[gene]
        pair = 100.0 * float(df.loc[gene, "frac_pairs_resolvable"])
        nearest = 100.0 * float(df.loc[gene, "frac_identifiable"])

        ax.plot([0, 100], [y, y], color=S.PALE_GREY, lw=0.55, zorder=0)
        ax.plot([0, pair], [y + dy, y + dy], color=color, lw=0.85, alpha=0.55, zorder=2)
        ax.plot(
            pair,
            y + dy,
            marker="o",
            ms=4.0,
            mfc="white",
            mec=color,
            mew=1.0,
            zorder=4,
        )
        ax.plot([0, nearest], [y - dy, y - dy], color=color, lw=1.0, zorder=2)
        ax.plot(
            nearest,
            y - dy,
            marker="D",
            ms=3.7,
            mfc=color,
            mec="white",
            mew=0.45,
            zorder=4,
        )

        for value, yy in ((pair, y + dy), (nearest, y - dy)):
            label = "0%" if value == 0 else f"{value:.1f}%"
            ax.text(
                value + 1.8,
                yy,
                label,
                ha="left",
                va="center",
                fontsize=5.7,
                color=S.INK,
                zorder=5,
            )

    labels = [f"{gene}\nn={int(df.loc[gene, 'n_var'])}" for gene in S.GENE_ORDER]
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, linespacing=1.05)
    S.gene_ticklabels(ax, "y")
    ax.set_xlim(-4, 103)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_ylim(-0.45, 3.45)
    ax.set_xlabel("Sibling-allele separation (%) · cosine-distance criterion", labelpad=2)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", pad=1.5)
    S.despine(ax, keep=("bottom",))

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            ls="-",
            lw=0.8,
            ms=4.0,
            color=S.GREY,
            mfc="white",
            mec=S.GREY,
            label="resolvable sibling pairs",
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            ls="-",
            lw=0.9,
            ms=3.6,
            color=S.GREY,
            mfc=S.GREY,
            mec=S.GREY,
            label="variants separable from nearest sibling",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.145, 0.99),
        ncol=2,
        columnspacing=1.25,
        handlelength=1.25,
        handletextpad=0.4,
        borderaxespad=0,
        fontsize=5.6,
    )
    S.save(fig, "fig3d_sibling_identification")


if __name__ == "__main__":
    main()
