"""Figure 3f (composite letter) - Effect size and sampling depth jointly determine the window.

Data-direct. Source: native-depth per-perturbation rows (fig3_data), plotting
variant-to-WT effect size (x, log) against the split-half ratio D_self/D_null (y),
with marker area proportional to effective cells per variant and colour by gene.
Rebuts the "TP53/KRAS are just low-depth" reading: they carry many cells yet sit
near R=1, whereas JAK1 has fewer cells but a wide (low-R) window. No regression is
fit; the panel shows the joint landscape only.

Layout notes (Nature Methods pass): the gene key is NOT repeated here. Genes are
labelled directly on their point clouds in the gene colour, with a hairline leader
where a cloud centre is crowded, so the only key in the panel is the marker-area
key for cells per variant (the one encoding not established anywhere else).

Run:  python fig3d_landscape.py  ->  fig3d_landscape.pdf (+ .png)
"""
from __future__ import annotations

import numpy as np
from matplotlib.lines import Line2D

import fig3_data as D
import nm_style as S

W_MM, H_MM = 60.3, 57.1

# Direct-label anchor (text position) and the cloud point the hairline leader runs
# to, in data units. Chosen once by eye from the plotted clouds; they annotate the
# data, they do not alter it.
LABEL_POS = {
    "KRAS":  ((0.64, 1.31), (1.00, 1.12)),
    "TP53":  ((0.64, 0.66), (1.32, 0.88)),
    "GATA1": ((6.4, 1.26), (3.4, 1.04)),
    "JAK1":  ((6.6, 0.36), (11.5, 0.13)),
}


def main() -> None:
    S.apply_rcparams()
    fig, ax = S.panel(W_MM, H_MM)

    for gene in S.GENE_ORDER:
        n = D.native_rankability(gene)
        es = n["effect_size"].to_numpy(float)
        ratio = n["ratio"].to_numpy(float)
        ncell = n["n_cells"].to_numpy(float)
        good = np.isfinite(es) & np.isfinite(ratio) & (es > 0)
        sizes = 3 + 22 * (ncell[good] - ncell[good].min()) / (np.ptp(ncell[good]) + 1e-9)
        # shape repeats the gene; marker area still encodes cells per variant
        ax.scatter(es[good], ratio[good], s=sizes, color=S.GENE_COLORS[gene],
                   marker=S.GENE_MARKERS[gene],
                   alpha=0.55, linewidths=0, zorder=3)

    ax.axhline(1.0, color=S.GREY, ls="--", lw=0.7, zorder=1)
    ax.set_xscale("log")
    # NOT "(energy D_null)": the plotted column is `effect_size`, and
    # effect_size / D_null ranges 1.01-6.65 across the four genes, so naming the
    # x-axis D_null told the reader it was the denominator of the y-axis and
    # manufactured a structural relationship that is not there.
    ax.set_xlabel("Variant-to-WT effect size")
    ax.set_ylabel(r"$D_\mathrm{self}/D_\mathrm{null}$")
    ax.set_xlim(0.62, 16.0)
    # plain tick labels: matplotlib's 10^n superscripts render at ~4.7 pt, below
    # the 5 pt house floor, and the decade range here is only ~1.5 decades
    ax.set_xticks([1, 2, 5, 10])
    ax.set_xticklabels(["1", "2", "5", "10"])
    ax.minorticks_off()
    ax.set_ylim(-0.03, 1.42)
    S.despine(ax)
    ax.tick_params(length=2.2)

    # ---- direct gene labels (replace the repeated gene legend) ----
    for gene, ((tx, ty), (lx, ly)) in LABEL_POS.items():
        col = S.GENE_COLORS[gene]
        ax.annotate(gene, xy=(lx, ly), xytext=(tx, ty), color=col,
                    fontsize=6, fontweight="bold", ha="left", va="center",
                    zorder=6,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color=col,
                                    shrinkA=2.0, shrinkB=2.0, alpha=0.8))

    # ---- the only key in this panel: marker area = cells per variant ----
    # neutral hexagon, deliberately not one of GENE_MARKERS: this key encodes cells
    # per variant, and reusing the TP53 circle would make one shape mean two things
    size_handles = [Line2D([0], [0], marker="h", ls="none", markeredgecolor="none",
                           markerfacecolor=S.GREY, markersize=ms, label=lab)
                    for ms, lab in [(1.9, "few"), (4.4, "many")]]
    leg = ax.legend(handles=size_handles, loc="lower left", title="cells per variant",
                    handletextpad=0.3, labelspacing=0.45, borderpad=0.0,
                    borderaxespad=0.2, title_fontsize=5.5, fontsize=5.5)
    leg.get_title().set_color(S.GREY)
    for t in leg.get_texts():
        t.set_color(S.GREY)

    S.save(fig, "fig3d_landscape")
    print("drawn: effect-size vs ratio, size=cells; genes", S.GENE_ORDER)


if __name__ == "__main__":
    main()
