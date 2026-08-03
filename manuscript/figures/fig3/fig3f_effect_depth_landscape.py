"""Fig. 3f — pseudobulk effect size, nominal depth and measurement window."""
from __future__ import annotations

import numpy as np
from matplotlib.lines import Line2D

import fig3_common as S

# Publication contract inherited from fig3_common: Arial/Helvetica sans-serif;
# svg.fonttype="none"; pdf.fonttype=42; exports .svg, .pdf, .png and .tiff at dpi=600.
# Final assembled figure target: width_mm = 183.
W_MM, H_MM = 78.0, 43.0

LABEL_POS = {
    "TP53": (0.72, 0.71),
    "KRAS": (0.72, 1.30),
    "GATA1": (5.8, 1.28),
    "JAK1": (6.1, 0.32),
}


def main() -> None:
    S.apply_style()
    df = S.nominal_depth_join()
    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0.15, 0.16, 0.81, 0.76])

    medians = {}
    for gene in S.GENE_ORDER:
        sub = df[df["dataset"] == gene]
        # Strictly positive effect-size guard is required for the log x-axis.
        good = (
            np.isfinite(sub["effect_size"])
            & np.isfinite(sub["window_ratio"])
            & np.isfinite(sub["n_cells_nominal"])
            & (sub["effect_size"] > 0)
        )
        use = sub.loc[good]
        ax.scatter(
            use["effect_size"],
            use["window_ratio"],
            s=S.size_from_nominal_cells(use["n_cells_nominal"]),
            color=S.GENE_COLORS[gene],
            alpha=0.32,
            linewidths=0,
            zorder=2,
        )
        med_x = float(use["effect_size"].median())
        med_y = float(use["window_ratio"].median())
        medians[gene] = (med_x, med_y)
        ax.plot(
            med_x,
            med_y,
            marker="D",
            ms=5.0,
            mfc=S.GENE_COLORS[gene],
            mec="white",
            mew=0.75,
            zorder=5,
        )

    ax.axhline(1.0, color=S.GREY, lw=0.7, ls=(0, (3, 2)), zorder=1)
    ax.text(15.7, 1.025, r"$R=1$", ha="right", va="bottom", fontsize=5.4, color=S.GREY)
    ax.set_xscale("log")
    ax.set_xlim(0.62, 16.5)
    ax.set_xticks([1, 2, 5, 10])
    ax.set_xticklabels(["1", "2", "5", "10"])
    ax.minorticks_off()
    ax.set_ylim(0, 1.42)
    ax.set_yticks(np.arange(0, 1.41, 0.2))
    ax.set_xlabel(r"Pseudobulk effect size, $\Vert\delta_v\Vert_2$", labelpad=2)
    ax.set_ylabel(r"Window ratio, $D_\mathrm{self}/D_\mathrm{null}$", labelpad=2)
    ax.tick_params(axis="x", pad=1.5)
    ax.tick_params(axis="y", pad=1.5)
    S.despine(ax)

    for gene, (tx, ty) in LABEL_POS.items():
        color = S.GENE_COLORS[gene]
        px, py = medians[gene]
        ax.annotate(
            gene,
            xy=(px, py),
            xytext=(tx, ty),
            color=color,
            fontsize=6.1,
            fontweight="bold",
            ha="left",
            va="center",
            arrowprops=dict(
                arrowstyle="-",
                lw=0.45,
                color=color,
                shrinkA=2.0,
                shrinkB=2.0,
                alpha=0.8,
            ),
            zorder=6,
        )

    size_values = [100, 1000, 10000]
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            ls="none",
            markeredgecolor="none",
            markerfacecolor=S.GREY,
            markersize=np.sqrt(S.size_from_nominal_cells([n])[0]),
            label=f"{n:,}",
        )
        for n in size_values
    ]
    legend = ax.legend(
        handles=handles,
        title="nominal cells per variant",
        loc="lower left",
        ncol=1,
        borderaxespad=0.2,
        handletextpad=0.5,
        labelspacing=0.35,
        fontsize=5.2,
        title_fontsize=5.2,
    )
    legend.get_title().set_color(S.GREY)
    for text in legend.get_texts():
        text.set_color(S.GREY)

    S.save(fig, "fig3f_effect_depth_landscape")


if __name__ == "__main__":
    main()
