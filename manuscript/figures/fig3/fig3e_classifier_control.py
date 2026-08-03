"""Fig. 3e — single-cell two-sample classifier control."""
from __future__ import annotations

from matplotlib.lines import Line2D

import fig3_common as S

# Publication contract inherited from fig3_common: Arial/Helvetica sans-serif;
# svg.fonttype="none"; pdf.fonttype=42; exports .svg, .pdf, .png and .tiff at dpi=600.
# Final assembled figure target: width_mm = 183.
W_MM, H_MM = 100.0, 43.0


def main() -> None:
    S.apply_style()
    df = S.read_csv("classifier_two_sample.csv").set_index("gene")
    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0.16, 0.24, 0.80, 0.60])

    y_positions = list(range(len(S.GENE_ORDER)))[::-1]
    dy = 0.12
    ax.axvline(0.5, color=S.GREY, lw=0.7, ls=(0, (3, 2)), zorder=1)
    for y, gene in zip(y_positions, S.GENE_ORDER):
        color = S.GENE_COLORS[gene]
        detect = float(df.loc[gene, "detect_med_auroc"])
        identify = float(df.loc[gene, "ident_med_auroc"])
        perm = float(df.loc[gene, "perm_med_auroc"])

        ax.plot([0.45, 1.0], [y, y], color=S.PALE_GREY, lw=0.55, zorder=0)
        ax.plot(
            [detect, identify],
            [y + dy, y - dy],
            color=color,
            lw=0.8,
            alpha=0.65,
            zorder=2,
        )
        ax.plot(
            detect,
            y + dy,
            marker="o",
            ms=4.2,
            mfc="white",
            mec=color,
            mew=1.0,
            zorder=4,
        )
        ax.plot(
            identify,
            y - dy,
            marker="D",
            ms=3.8,
            mfc=color,
            mec="white",
            mew=0.45,
            zorder=4,
        )
        ax.plot(
            perm,
            y,
            marker="x",
            ms=3.5,
            color=S.GREY,
            mew=0.75,
            zorder=5,
        )

        if gene == "JAK1":
            ax.text(
                detect - 0.012,
                y + dy + 0.13,
                f"{detect:.2f}",
                ha="center",
                va="bottom",
                fontsize=5.7,
                color=color,
            )
            ax.text(
                identify - 0.004,
                y - dy - 0.13,
                f"{identify:.2f}",
                ha="center",
                va="top",
                fontsize=5.7,
                color=color,
            )

    labels = [f"{gene}\nn={int(df.loc[gene, 'n_var'])}" for gene in S.GENE_ORDER]
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, linespacing=1.05)
    S.gene_ticklabels(ax, "y")
    ax.set_xlim(0.45, 1.0)
    ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_ylim(-0.45, 3.45)
    ax.set_xlabel("Cross-validated two-sample AUROC", labelpad=2)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", pad=1.5)
    S.despine(ax, keep=("bottom",))

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            ls="none",
            ms=4.1,
            mfc="white",
            mec=S.GREY,
            mew=1.0,
            label="detection: variant vs WT",
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            ls="none",
            ms=3.7,
            mfc=S.GREY,
            mec=S.GREY,
            label="identification: sibling vs sibling",
        ),
        Line2D(
            [0],
            [0],
            marker="x",
            ls="none",
            ms=3.4,
            color=S.GREY,
            label="permutation median",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.145, 0.99),
        ncol=3,
        columnspacing=0.85,
        handlelength=0.65,
        handletextpad=0.3,
        borderaxespad=0,
        fontsize=5.25,
    )
    ax.text(
        0.505,
        -0.38,
        "chance",
        ha="left",
        va="bottom",
        fontsize=5.2,
        color=S.GREY,
    )
    S.save(fig, "fig3e_classifier_control")


if __name__ == "__main__":
    main()
