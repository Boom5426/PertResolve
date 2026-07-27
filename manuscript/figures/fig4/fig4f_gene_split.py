"""Figure 4f - Direction/magnitude dissociation holds across gene x split.

Data-direct. Source (committed): results/results_v4_exttheta.csv (per-variant grid).
Aggregation: for each (gene, split) take the mean pearson_delta and mean PDS_cos
over the 18 in-house feature-model heads (Ridge/Lasso/RF/GBoost/KNN/MLP x
theta/esm/esm+theta), one scatter point per gene-split combination.

Verified firsthand (17 gene-split combinations):
  pearson_delta > PDS_cos in 15 of 17 combos (points above the y=x diagonal)
  the 2 exceptions are JAK1 Mechanistic (split3) and KRAS Mechanistic (split3)
  GATA1 Low-depth (split5): PDS_cos = 0.186 (deepest dissociation; labelled)

Message: direction can be recovered (high pearson_delta) while per-variant
identity is not (PDS_cos near chance); this dissociation holds across nearly all
gene-split combinations, not a single split or gene.

Legend economy (Nature Methods pass): the detached four-entry gene key is replaced
by direct labels placed next to each gene's points, so the figure carries no
repeated gene legend. Gene hues are the house GENE_COLORS and appear only in this
panel of Fig. 4; every other Fig. 4 panel is drawn in house neutrals, so no hue is
doing two jobs. Axis names match the metric names used in Fig. 4a,b,e.

Run:  python fig4f_gene_split.py  ->  fig4f_gene_split.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

# direct-label anchors, in data coordinates, chosen to sit beside each gene's
# points without overlapping any marker or another label.
GENE_LABEL_XY = {
    "TP53":  (0.42, 0.795, "center", "bottom"),
    "KRAS":  (0.65, 0.742, "center", "bottom"),
    "JAK1":  (0.595, 0.665, "left", "center"),
    "GATA1": (0.595, 0.575, "left", "center"),
}


def main() -> None:
    S.apply_rcparams()
    df = D.exttheta()
    inh = df[df.method.apply(D.is_inhouse)].copy()
    g = (inh.groupby(["gene", "split"])
             .agg(pds=("PDS_cos", "mean"), pear=("pearson_delta", "mean"))
             .reset_index())

    n = len(g)
    n_above = int((g.pear > g.pds).sum())

    fig, ax = S.panel(52, 47.3)

    # y=x diagonal (dissociation boundary): points above => pearson > PDS.
    # The line is left unlabelled: the in-panel note and the caption already say
    # what crossing it means, and a rotated two-character label is the smallest
    # type on the page.
    lim_lo, lim_hi = 0.0, 0.82
    ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], ls="--", lw=0.6,
            color=S.LIGHT_GREY, zorder=1)

    # chance reference for PDS (0.5) as a light vertical guide
    ax.axvline(0.5, color=S.LIGHT_GREY, lw=0.5, ls=":", zorder=0)
    ax.text(0.507, 0.045, "chance", fontsize=5.0,
            color=S.GREY, ha="left", va="bottom", rotation=90)

    for gene in S.GENE_ORDER:
        sub = g[g.gene == gene]
        ax.scatter(sub.pds, sub.pear, s=16, c=S.GENE_COLORS[gene],
                   edgecolors="white", linewidths=0.4, zorder=3)
        x, y, ha, va = GENE_LABEL_XY[gene]
        ax.text(x, y, gene, fontsize=6.0, color=S.GENE_COLORS[gene],
                ha=ha, va=va, zorder=4)

    # label the deepest dissociation: GATA1 Low-depth (split5)
    gp = g[(g.gene == "GATA1") & (g.split == "split5")].iloc[0]
    ax.annotate("GATA1 Low-depth\nPDS 0.19", xy=(gp.pds, gp.pear),
                xytext=(gp.pds + 0.030, gp.pear - 0.185), fontsize=5.0,
                color=S.INK, ha="left", va="center", linespacing=1.2,
                arrowprops=dict(arrowstyle="-", lw=0.5, color=S.GREY,
                                shrinkA=0, shrinkB=2.5))

    ax.set_xlim(lim_lo, lim_hi)
    ax.set_ylim(lim_lo, lim_hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("PDS (cosine)")
    ax.set_ylabel(r"Pearson-$\delta$ (direction)")
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8])
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
    S.despine(ax)
    ax.tick_params(length=2.2)

    ax.text(0.02, 0.985,
            f"Pearson-$\\delta$ > PDS\nin {n_above} of {n}\ngene-split combos",
            transform=ax.transAxes, fontsize=5.2, color=S.INK,
            va="top", ha="left", linespacing=1.35)

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig4f_gene_split"))
    print(f"n={n} combos  pearson>PDS in {n_above}/{n}  "
          f"GATA1 split5 PDS={gp.pds:.3f} pear={gp.pear:.3f}")


if __name__ == "__main__":
    main()
