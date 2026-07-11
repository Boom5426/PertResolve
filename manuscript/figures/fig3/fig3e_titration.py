"""Figure 3e - Depth titration reveals gene-dependent gains.

Data-direct. Source: results/split_half_power_curve.csv, PCA-50 space (canonical).
Fraction of variants detectable vs cells per half (50/100/150/300). JAK1 is already
saturated at 50 cells/half; TP53, KRAS and GATA1 improve only modestly and plateau
below 75%. (The committed curve stores point estimates only, no bootstrap ribbon;
TP53 has no 300-cell rung.)
Message: adding cells helps only when the allele-specific effect window is wide enough.

Run:  python fig3e_titration.py  ->  fig3e_titration.pdf (+ .png)
"""
from __future__ import annotations

import fig3_data as D
import nm_style as S


def main() -> None:
    S.apply_rcparams()
    pw = D.power_curve("pca50")

    fig, ax = S.panel(56, 46)
    for gene in S.GENE_ORDER:
        s = pw[pw.gene == gene].sort_values("n_sub")
        ax.plot(s.n_sub, s.frac_detectable, "-o", color=S.GENE_COLORS[gene],
                lw=1.0, ms=3.0, label=gene, zorder=3)

    ax.axhline(0.75, color=S.GREY, ls=":", lw=0.6, zorder=1)
    ax.text(305, 0.75, "75%", fontsize=5, color=S.GREY, va="center", ha="left")
    ax.set_xlabel("Cells per half")
    ax.set_ylabel("Fraction detectable")
    ax.set_xticks([50, 100, 150, 300])
    ax.set_ylim(0.0, 1.05)
    ax.set_xlim(35, 330)
    S.despine(ax)
    ax.tick_params(length=2.2)
    leg = ax.legend(loc="center right", handletextpad=0.4, labelspacing=0.3,
                    borderpad=0.2)
    for t, g in zip(leg.get_texts(), S.GENE_ORDER):
        t.set_color(S.GENE_COLORS[g])

    S.save(fig, "fig3e_titration")
    print(pw.pivot_table(index="gene", columns="n_sub", values="frac_detectable")
          .reindex(S.GENE_ORDER).round(3).to_string())


if __name__ == "__main__":
    main()
