"""Figure 3e (composite letter) - Model-free two-sample classifier AUROC.

Data-direct. Source: results/canonical/classifier_two_sample.csv via D.classifier().
Per gene, a supervised single-cell logistic classifier is trained/evaluated two ways:
  - detection : variant cells vs wild-type cells (detect_med_auroc, open marker)
  - identification : one sibling allele vs another sibling allele (ident_med_auroc, filled)
A permutation control (perm_med_auroc ~0.50) and the chance line (0.50) bound the null.

Message: a supervised classifier separates variant-from-WT and sibling-from-sibling only
for JAK1; for TP53/KRAS/GATA1 both stay at chance, so detection (vs WT) is not the same
as identification (which allele).

Numbers plotted (median AUROC across variants / sibling pairs), read from the CSV:
  gene   detect  ident   perm
  TP53   0.533   0.493   0.501
  KRAS   0.498   0.498   0.501
  GATA1  0.508   0.520   0.500
  JAK1   0.965   0.873   0.503

Layout notes (Nature Methods pass): drawn at its final placed size so on-page type is
~6 pt (composite scale ~1.0) and it shares a common top edge and baseline with panel d;
the open/filled marker key sits on one quiet header line, and gene identity is carried
by the coloured y-axis labels (no gene legend is repeated here).

Run:  python fig3h_classifier.py  ->  fig3h_classifier.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

W_MM, H_MM = 96.0, 31.9


def main() -> None:
    S.apply_rcparams()
    cls = D.classifier().set_index("gene")

    fig, ax = S.panel(W_MM, H_MM)
    genes = S.GENE_ORDER
    ys = list(range(len(genes)))

    # chance line at 0.50
    ax.axvline(0.50, color=S.GREY, lw=0.6, ls=(0, (3, 2)), zorder=1)

    for y, gene in zip(ys, genes):
        det = float(cls.loc[gene, "detect_med_auroc"])
        idn = float(cls.loc[gene, "ident_med_auroc"])
        color = S.GENE_COLORS[gene]
        lo, hi = sorted((det, idn))
        # connector
        ax.plot([lo, hi], [y, y], color=color, lw=1.1, zorder=2,
                solid_capstyle="round")
        # detection = open marker (variant vs WT)
        ax.plot(det, y, marker="o", ms=4.2, mfc="white", mec=color, mew=1.0,
                zorder=4, clip_on=False)
        # identification = filled marker (sibling vs sibling)
        ax.plot(idn, y, marker="o", ms=4.2, mfc=color, mec=color, mew=1.0,
                zorder=4, clip_on=False)

    ax.set_yticks(ys)
    ax.set_yticklabels(genes)
    for t, g in zip(ax.get_yticklabels(), genes):
        t.set_color(S.GENE_COLORS[g])
        t.set_fontweight("bold")
    ax.invert_yaxis()

    ax.set_xlim(0.45, 1.0)
    ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_xlabel("Two-sample classifier AUROC")

    S.despine(ax, keep=("bottom",))
    ax.tick_params(length=2.2)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(genes) - 0.42, -0.88)

    # chance annotation, direct on the line inside the plot
    ax.text(0.507, 3.30, "chance", va="center", ha="left", fontsize=5.5,
            color=S.GREY)

    # one quiet marker key on a single header line (neutral grey: it names the
    # marker STYLE, not the gene colour, which the y-axis labels already carry)
    handles = [
        Line2D([0], [0], marker="o", ls="none", ms=4.2, mfc="white",
               mec=S.GREY, mew=1.0, label="detection (vs wild type)"),
        Line2D([0], [0], marker="o", ls="none", ms=4.2, mfc=S.GREY,
               mec=S.GREY, mew=1.0, label="identification (sibling vs sibling)"),
    ]
    leg = ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.005),
                    ncol=2, handlelength=1.0, columnspacing=1.4,
                    handletextpad=0.35, borderaxespad=0.0, borderpad=0.0,
                    fontsize=5.5)
    for t in leg.get_texts():
        t.set_color(S.GREY)
    leg.set_zorder(6)

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig3h_classifier"))

    for g in genes:
        print(f"{g:6} detect={cls.loc[g,'detect_med_auroc']:.3f}  "
              f"ident={cls.loc[g,'ident_med_auroc']:.3f}  "
              f"perm={cls.loc[g,'perm_med_auroc']:.3f}")


if __name__ == "__main__":
    main()
