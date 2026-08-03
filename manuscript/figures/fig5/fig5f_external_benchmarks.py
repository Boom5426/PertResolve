"""Fig. 5f — seven-level rank recovery across nine matched-depth datasets."""
from __future__ import annotations

W_MM, H_MM = 68.0, 34.0

from matplotlib.lines import Line2D

import fig5_common as S

ALLELE = {
    "allele_TP53": "TP53",
    "allele_KRAS": "KRAS",
    "allele_GATA1": "GATA1",
    "allele_JAK1": "JAK1",
}
DISPLAY = {
    "sciPlex": "sci-Plex",
    "Norman": "Norman",
    "Replogle": "Replogle",
    "Adamson": "Adamson",
    "VCC": "VCC",
}
# dx, dy, horizontal alignment
OFFSETS = {
    "allele_TP53": (0.013, 0.060, "left"),
    "allele_KRAS": (0.013, -0.018, "left"),
    "allele_GATA1": (-0.012, 0.058, "right"),
    "allele_JAK1": (0.010, -0.055, "left"),
    "sciPlex": (0.000, -0.075, "center"),
    "VCC": (-0.010, -0.060, "right"),
    "Replogle": (-0.005, 0.060, "right"),
    "Adamson": (-0.010, -0.055, "right"),
    "Norman": (0.008, 0.060, "left"),
}


def main() -> None:
    S.apply_style()
    data = S.read_csv("benchmark_resolution_summary.csv")

    fig, ax = S.panel(W_MM, H_MM)
    fig.subplots_adjust(left=0.235, right=0.965, bottom=0.29, top=0.95)

    ax.axvline(0.5, color=S.MID_GREY, lw=0.65, ls=(0, (3, 2)), zorder=0)
    ax.axhspan(0.90, 1.055, color=S.PALE_GREY, zorder=0)
    ax.axhspan(-0.04, 0.50, color="#F7F7F7", zorder=0)
    ax.text(0.468, 0.845, "reliable full-order recovery", fontsize=5.2,
            color=S.GREY, va="center")
    ax.text(0.700, 0.105, "ranking unstable", fontsize=5.2,
            color=S.GREY, ha="left", va="center")

    for _, row in data.iterrows():
        dataset = row["dataset"]
        x = float(row["oracle_ceiling"])
        y = float(row["resolution_P_recover_order"])
        display_y = {
            "allele_TP53": 0.050,
            "allele_KRAS": 0.005,
        }.get(dataset, y)
        n = float(row["n_perturbations"])
        if dataset in ALLELE:
            gene = ALLELE[dataset]
            color = S.GENE_COLORS[gene]
            label = gene
        else:
            color = S.GREY
            label = DISPLAY[dataset]
        if display_y != y:
            ax.plot([x, x], [y, display_y], color=color, lw=0.65,
                    alpha=0.7, zorder=2)
            ax.scatter([x], [y], s=5.0, color=color, edgecolor="none", zorder=3)
        ax.scatter([x], [display_y], s=S.point_size(n), color=color,
                   edgecolor="white", linewidths=0.55, zorder=4)
        dx, dy, ha = OFFSETS[dataset]
        ax.text(x + dx, display_y + dy, label, fontsize=5.6,
                color=S.INK if dataset not in ALLELE else color,
                ha=ha, va="center", zorder=5)

    size_values = [25, 100, 400]
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=S.MID_GREY,
               markeredgecolor="white", markeredgewidth=0.4,
               markersize=(S.point_size(value) ** 0.5) * 0.86,
               label=str(value))
        for value in size_values
    ]
    # Inside the axes, lower right: that corner of the "ranking unstable" band
    # carries no markers, so the key costs no plot area and frees the outside
    # band the panel used to reserve for it.
    legend = ax.legend(handles=handles, title="perturbations",
                       loc="lower right", borderaxespad=0.35,
                       ncol=1, borderpad=0.25, labelspacing=0.38,
                       handletextpad=0.5)
    legend.get_title().set_fontsize(5.5)
    legend.get_title().set_color(S.GREY)

    ax.set_xlim(0.46, 0.98)
    ax.set_ylim(-0.04, 1.055)
    ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9])
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    # "matched 50 cells per half" is in the caption; at 68 mm the full label
    # is 58 mm of text on a 47 mm axis and runs off the panel.
    ax.set_xlabel("Replicate ceiling (PDS)")
    # Shortened for a 21 mm axis: the full wording rotated to 25 mm and overran
    # the 32 mm canvas. The caption gives it in full.
    ax.set_ylabel("P(recover order)")
    S.despine(ax)

    S.save(fig, "fig5f_external_benchmarks")


if __name__ == "__main__":
    main()
