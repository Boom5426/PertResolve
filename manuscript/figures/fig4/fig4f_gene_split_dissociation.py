"""Fig. 4f — direction–ranking dissociation across gene-by-split combinations."""
from __future__ import annotations

from matplotlib.lines import Line2D

import fig4_common as S

# Static QA contract: Arial/Helvetica sans-serif; svg.fonttype="none";
# pdf.fonttype=42; exports .svg, .pdf, .png and .tiff at dpi=600.
# Final assembled figure target: width_mm = 183.
W_MM, H_MM = 89.0, 45.0

SPLIT_MARKERS = {
    "split1": "o",
    "split2": "s",
    "split3": "^",
    "split5": "D",
    "split6": "P",
}


def main() -> None:
    S.apply_style()
    df = S.breakdown()
    inhouse = df[df["method"].map(S.is_inhouse)]
    summary = (
        inhouse.groupby(["gene", "split"])
        .agg(pds=("PDS_cos", "mean"), direction=("pearson_delta", "mean"))
        .reset_index()
    )
    n_above = int((summary["direction"] > summary["pds"]).sum())
    n = len(summary)

    fig = S.figure(W_MM, H_MM)
    # Plot on the left, keys in a column on the right. Stacked above the
    # plot the two key rows cost 12 mm of the 45 mm panel and forced the
    # square axes down to a third of the width; beside it they cost none.
    ax = fig.add_axes([0.125, 0.205, 0.50, 0.755])
    ax.set_anchor("W")
    ax.plot([0, 0.82], [0, 0.82], color=S.LIGHT_GREY, lw=0.7, ls=(0, (3, 2)))
    ax.axvline(0.50, color=S.LIGHT_GREY, lw=0.6, ls=(0, (1.5, 2)))
    for _, row in summary.iterrows():
        ax.plot(
            row["pds"],
            row["direction"],
            marker=SPLIT_MARKERS[row["split"]],
            ms=4.2,
            mfc=S.GENE_COLORS[row["gene"]],
            mec="white",
            mew=0.5,
            ls="none",
            zorder=4,
        )

    exceptions = summary[summary["direction"] <= summary["pds"]]
    ax.scatter(
        exceptions["pds"],
        exceptions["direction"],
        s=38,
        facecolors="none",
        edgecolors=S.AMBER,
        linewidths=0.9,
        zorder=5,
    )

    low_depth = summary[
        (summary["gene"] == "GATA1") & (summary["split"] == "split5")
    ].iloc[0]
    ax.annotate(
        "GATA1 Low-depth\nPDS 0.19",
        xy=(low_depth["pds"], low_depth["direction"]),
        xytext=(0.06, 0.37),
        ha="left",
        va="center",
        fontsize=5.0,
        color=S.INK,
        arrowprops=dict(arrowstyle="-", lw=0.5, color=S.GREY),
    )
    ax.text(
        0.507,
        0.04,
        "PDS chance",
        rotation=90,
        ha="left",
        va="bottom",
        fontsize=5.0,
        color=S.GREY,
    )

    ax.set_xlim(0, 0.82)
    ax.set_ylim(0, 0.82)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.set_xlabel("PDS (cosine)")
    ax.set_ylabel(r"Pearson-$\delta$ (direction)")
    ax.tick_params(axis="both", pad=1.2)
    S.despine(ax)

    gene_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            ls="none",
            ms=3.8,
            mfc=S.GENE_COLORS[gene],
            mec="white",
            mew=0.4,
            label=gene,
        )
        for gene in S.GENE_ORDER
    ]
    split_handles = [
        Line2D(
            [0],
            [0],
            marker=SPLIT_MARKERS[split],
            ls="none",
            ms=3.5,
            mfc=S.MID_GREY,
            mec="white",
            mew=0.4,
            label={
                "split1": "Random",
                "split2": "Positional",
                "split3": "Mechanistic",
                "split5": "Low-depth",
                "split6": "Compat.",
            }[split],
        )
        for split in S.SPLIT_ORDER
    ]
    fig.legend(
        handles=gene_handles,
        loc="upper left",
        bbox_to_anchor=(0.615, 0.845),
        ncol=1,
        labelspacing=0.42,
        handletextpad=0.3,
        borderaxespad=0,
        borderpad=0,
        fontsize=5.0,
    )
    fig.legend(
        handles=split_handles,
        loc="upper left",
        bbox_to_anchor=(0.805, 0.845),
        ncol=1,
        labelspacing=0.42,
        handletextpad=0.3,
        borderaxespad=0,
        borderpad=0,
        fontsize=5.0,
    )
    fig.text(
        0.615,
        0.42,
        "amber outline:\n" r"Pearson-$\delta$ $\leq$ PDS",
        ha="left",
        va="top",
        fontsize=5.0,
        color=S.GREY,
        linespacing=1.25,
    )
    S.save(fig, "fig4f_gene_split_dissociation")


if __name__ == "__main__":
    main()
