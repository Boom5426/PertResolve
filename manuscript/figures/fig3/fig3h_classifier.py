"""Figure 3h - Model-free two-sample classifier AUROC (detection vs identification).

Data-direct. Source: results/_remote/unified/classifier_two_sample.csv via D.classifier().
Per gene, a supervised single-cell logistic classifier is trained/evaluated two ways:
  - detection : variant cells vs wild-type cells (detect_med_auroc, open marker)
  - identification : one sibling allele vs another sibling allele (ident_med_auroc, filled)
A permutation control (perm_med_auroc ~0.50) and the chance line (0.50) bound the null.

Message: a supervised classifier separates variant-from-WT and sibling-from-sibling only
for JAK1; for TP53/KRAS/GATA1 both stay at chance, so detection (vs WT) is not the same
as identification (which allele).

Numbers plotted (median AUROC across variants / sibling pairs):
  gene   detect  ident   perm
  TP53   0.533   0.493   0.501
  KRAS   0.494   0.496   0.504
  GATA1  0.508   0.520   0.500
  JAK1   0.965   0.873   0.503

Run:  python fig3h_classifier.py  ->  fig3h_classifier.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D


def main() -> None:
    S.apply_rcparams()
    cls = D.classifier().set_index("gene")

    fig, ax = S.panel(56, 44)
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
        ax.plot(det, y, marker="o", ms=5.0, mfc="white", mec=color, mew=1.1,
                zorder=4, clip_on=False)
        # identification = filled marker (sibling vs sibling)
        ax.plot(idn, y, marker="o", ms=5.0, mfc=color, mec=color, mew=1.1,
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

    # chance annotation
    ax.text(0.50, -0.62, "chance", va="bottom", ha="center", fontsize=5.5,
            color=S.GREY)

    S.despine(ax, keep=("bottom",))
    ax.tick_params(length=2.2)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(genes) - 0.4, -0.9)

    # legend: open = detection (vs WT), filled = identification (sibling)
    lx = 0.615
    ly0, dy = 0.30, 0.55
    ax.plot(lx, ly0, marker="o", ms=5.0, mfc="white", mec=S.INK, mew=1.1,
            clip_on=False, zorder=5)
    ax.text(lx + 0.022, ly0, "detection (vs WT)", va="center", ha="left",
            fontsize=5.8, color=S.INK)
    ax.plot(lx, ly0 + dy, marker="o", ms=5.0, mfc=S.INK, mec=S.INK, mew=1.1,
            clip_on=False, zorder=5)
    ax.text(lx + 0.022, ly0 + dy, "identification (sibling vs sibling)",
            va="center", ha="left", fontsize=5.8, color=S.INK)

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig3h_classifier"))

    for g in genes:
        print(f"{g:6} detect={cls.loc[g,'detect_med_auroc']:.3f}  "
              f"ident={cls.loc[g,'ident_med_auroc']:.3f}  "
              f"perm={cls.loc[g,'perm_med_auroc']:.3f}")


if __name__ == "__main__":
    main()
