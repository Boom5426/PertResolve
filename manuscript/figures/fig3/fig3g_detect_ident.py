"""Figure 3g - Detection vs identification among sibling alleles (four genes).

Data-direct. Source: results/pairwise_resolvability.csv (D.pairwise()).
Two paired horizontal bars per gene, expressed as %:
  - "sibling pairs resolvable"  = frac_pairs_resolvable
  - "variants identifiable"     = frac_identifiable

Committed values:
  gene    frac_pairs_resolvable   frac_identifiable   median_nn_dist  median_Dself
  TP53    0.000                   0.000               0.230           0.965
  KRAS    0.000                   0.000               0.376           0.979
  GATA1   0.038                   0.000               0.461           0.981
  JAK1    0.785                   0.154               0.197           0.275

Message: allele-allele identification is ~0 for TP53/KRAS/GATA1 and only high
for JAK1 (78.5% of sibling pairs resolvable, 15.4% of variants identifiable),
because JAK1 sibling separation (nn_dist) is small relative to its self-noise
floor (D_self ~0.28) whereas the other three genes sit near D_self ~0.97.

Run:  python fig3g_detect_ident.py  ->  fig3g_detect_ident.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D


def main() -> None:
    S.apply_rcparams()
    pw = D.pairwise().set_index("gene")

    fig, ax = S.panel(56, 44)
    genes = S.GENE_ORDER

    # two bars per gene: resolvable (solid) and identifiable (hatched/lighter)
    group_h = 0.74          # total vertical span of a gene's pair
    bar_h = group_h / 2.0
    centers = list(range(len(genes)))

    for c, gene in zip(centers, genes):
        color = S.GENE_COLORS[gene]
        resolv = float(pw.loc[gene, "frac_pairs_resolvable"]) * 100.0
        ident = float(pw.loc[gene, "frac_identifiable"]) * 100.0

        y_res = c - bar_h / 2.0
        y_id = c + bar_h / 2.0

        # resolvable: filled bar
        ax.barh(y_res, resolv, height=bar_h * 0.9, color=color, alpha=0.9,
                edgecolor="none", zorder=3)
        # identifiable: hollow bar with the gene colour edge + light fill
        ax.barh(y_id, ident, height=bar_h * 0.9, facecolor="white",
                edgecolor=color, linewidth=0.7, hatch="////", zorder=3)

        # value labels just past bar end (or just past axis origin when 0)
        for y, val in ((y_res, resolv), (y_id, ident)):
            lab = f"{val:.1f}%" if val not in (0.0, 100.0) else "0%"
            ax.text(max(val, 0.0) + 1.5, y, lab, va="center", ha="left",
                    fontsize=5.5, color=S.INK)

    ax.set_yticks(centers)
    ax.set_yticklabels(genes)
    for t, g in zip(ax.get_yticklabels(), genes):
        t.set_color(S.GENE_COLORS[g])
        t.set_fontweight("bold")
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Fraction (%)")
    S.despine(ax, keep=("bottom",))
    ax.tick_params(length=2.2)
    ax.tick_params(axis="y", length=0)

    # legend as neutral swatches (grey) so it explains bar STYLE not gene colour
    import matplotlib.patches as mpatches
    solid = mpatches.Patch(facecolor=S.GREY, edgecolor="none",
                           label="sibling pairs resolvable")
    hollow = mpatches.Patch(facecolor="white", edgecolor=S.GREY, linewidth=0.7,
                            hatch="////", label="variants identifiable")
    leg = ax.legend(handles=[solid, hollow], loc="upper right",
                    bbox_to_anchor=(1.0, 1.0), handlelength=1.1,
                    handleheight=1.0, borderaxespad=0.0, labelspacing=0.5,
                    fontsize=5.5)
    leg.set_zorder(5)

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig3g_detect_ident"))

    # ---- report numbers ----
    for g in genes:
        print(f"{g:6} resolvable={pw.loc[g,'frac_pairs_resolvable']*100:5.1f}%  "
              f"identifiable={pw.loc[g,'frac_identifiable']*100:5.1f}%  "
              f"nn_dist={pw.loc[g,'median_nn_dist']:.3f}  "
              f"D_self={pw.loc[g,'median_Dself']:.3f}")


if __name__ == "__main__":
    main()
