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

Run:  python fig4f_gene_split.py  ->  fig4f_gene_split.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D


def main() -> None:
    S.apply_rcparams()
    df = D.exttheta()
    inh = df[df.method.apply(D.is_inhouse)].copy()
    g = (inh.groupby(["gene", "split"])
             .agg(pds=("PDS_cos", "mean"), pear=("pearson_delta", "mean"))
             .reset_index())

    n = len(g)
    n_above = int((g.pear > g.pds).sum())

    fig, ax = S.panel(56, 46)

    # y=x diagonal (dissociation boundary): points above => pearson > PDS
    lim_lo, lim_hi = 0.0, 0.82
    ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], ls="--", lw=0.6,
            color=S.GREY, zorder=1)
    ax.text(0.60, 0.615, "y = x", fontsize=5, color=S.GREY,
            ha="left", va="center", rotation=45, rotation_mode="anchor")

    # chance reference for PDS (0.5) as a light vertical guide
    ax.axvline(0.5, color=S.LIGHT_GREY, lw=0.5, ls=":", zorder=0)
    ax.text(0.505, 0.055, "PDS = 0.5 (chance)", fontsize=4.6,
            color=S.GREY, ha="left", va="bottom", rotation=90)

    for gene in S.GENE_ORDER:
        sub = g[g.gene == gene]
        ax.scatter(sub.pds, sub.pear, s=16, c=S.GENE_COLORS[gene],
                   edgecolors="white", linewidths=0.4, zorder=3, label=gene)

    # label the deepest dissociation: GATA1 Low-depth (split5)
    gp = g[(g.gene == "GATA1") & (g.split == "split5")].iloc[0]
    ax.annotate("GATA1\nLow-depth\nPDS 0.19", xy=(gp.pds, gp.pear),
                xytext=(gp.pds + 0.10, gp.pear - 0.145), fontsize=4.8,
                color=S.INK, ha="left", va="center", linespacing=1.15,
                arrowprops=dict(arrowstyle="-", lw=0.5, color=S.GREY,
                                shrinkA=0, shrinkB=2.5))

    ax.set_xlim(lim_lo, lim_hi)
    ax.set_ylim(lim_lo, lim_hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("PDS (cosine): per-variant identity")
    ax.set_ylabel("pearson_delta: direction")
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8])
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
    S.despine(ax)
    ax.tick_params(length=2.2)

    leg = ax.legend(loc="lower right", handletextpad=0.3, borderpad=0.3,
                    labelspacing=0.25, fontsize=6, markerscale=1.0)
    for h in leg.legend_handles:
        h.set_edgecolor("white")

    ax.text(0.03, 0.975,
            f"pearson > PDS in\n{n_above} of {n} gene-split combos",
            transform=ax.transAxes, fontsize=5.2, color=S.INK,
            va="top", ha="left", linespacing=1.35)

    S.save(fig, "fig4f_gene_split")
    print(f"n={n} combos  pearson>PDS in {n_above}/{n}  "
          f"GATA1 split5 PDS={gp.pds:.3f} pear={gp.pear:.3f}")


if __name__ == "__main__":
    main()
