"""Figure 3f - Fraction of variants rankable at native depth (four allele genes).

Data-direct. Source: results/unrankable_canonical.json (rankable = 100 - un-rankable),
cross-checked against fig3_data.rankable_fraction. Values (95% CI):
  TP53 0% (0-0), KRAS 0% (0-0), GATA1 2.4% (0.8-4.3), JAK1 90% (75-100).
This is the allele-level summary only; the gene-level atlases (Replogle/Norman/
Adamson) belong in Fig 6 / Extended Data, not here.

Run:  python fig3f_rankable.py  ->  fig3f_rankable.pdf (+ .png)
"""
from __future__ import annotations

import fig3_data as D
import nm_style as S


def main() -> None:
    S.apply_rcparams()
    uc = D.unrankable_canonical()

    fig, ax = S.panel(52, 38)
    genes = S.GENE_ORDER
    ys = range(len(genes))
    for y, gene in zip(ys, genes):
        rank = 100.0 - uc[gene]["ur"]
        lo = 100.0 - uc[gene]["ci_hi"]
        hi = 100.0 - uc[gene]["ci_lo"]
        color = S.GENE_COLORS[gene]
        ax.barh(y, rank, height=0.6, color=color, alpha=0.85, edgecolor="none", zorder=3)
        ax.plot([lo, hi], [y, y], color=S.INK, lw=0.9, zorder=4, solid_capstyle="round")
        lab = f"{rank:.1f}%" if rank not in (0.0, 100.0) else f"{rank:.0f}%"
        ax.text(104, y, lab, va="center", ha="left", fontsize=6, color=S.INK)  # right column

    ax.set_yticks(list(ys))
    ax.set_yticklabels(genes)
    for t, g in zip(ax.get_yticklabels(), genes):
        t.set_color(S.GENE_COLORS[g]); t.set_fontweight("bold")
    ax.invert_yaxis()
    ax.set_xlim(0, 122)
    ax.set_xlabel("Variants rankable at native depth (%)")
    ax.set_xticks([0, 25, 50, 75, 100])
    S.despine(ax, keep=("bottom",))
    ax.tick_params(length=2.2)
    ax.tick_params(axis="y", length=0)

    S.save(fig, "fig3f_rankable")
    for g in genes:
        print(f"{g:6} rankable={100-uc[g]['ur']:.1f}%  n={uc[g]['n']}")


if __name__ == "__main__":
    main()
