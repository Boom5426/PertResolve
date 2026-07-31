"""Figure 3g (composite letter) - Depth titration reveals gene-dependent gains.

Data-direct. Source: results/split_half_power_curve.csv, PCA-50 space (canonical).
Fraction of variants detectable vs cells per half (50/100/150/300). JAK1 is already
saturated at 50 cells/half; TP53, KRAS and GATA1 improve only modestly and plateau
below 75%. (The committed curve stores point estimates only, no bootstrap ribbon;
TP53 has no 300-cell rung.)
Message: adding cells helps only when the allele-specific effect window is wide enough.

Layout notes (Nature Methods pass): the gene key is NOT repeated here. Each curve is
labelled directly at its end in the gene colour, which also makes the panel readable
in greyscale (position identifies the curve, colour only reinforces it).

Run:  python fig3e_titration.py  ->  fig3e_titration.pdf (+ .png)
"""
from __future__ import annotations

import fig3_data as D
import nm_style as S

W_MM, H_MM = 54.6, 58.5


def main() -> None:
    S.apply_rcparams()
    pw = D.power_curve("pca50")

    fig, ax = S.panel(W_MM, H_MM)
    ends = {}
    for gene in S.GENE_ORDER:
        s = pw[pw.gene == gene].sort_values("n_sub")
        ax.plot(s.n_sub, s.frac_detectable, "-o", color=S.GENE_COLORS[gene],
                marker=S.GENE_MARKERS[gene], lw=1.0, ms=2.6, zorder=3)
        # n at the deepest rung is carried to the end label: the number of variants
        # that still HAVE that rung falls steeply with depth (JAK1 reaches 300 cells
        # per half with a single variant), and a curve endpoint resting on n = 1
        # must not read like the endpoints beside it.
        ends[gene] = (float(s.n_sub.iloc[-1]), float(s.frac_detectable.iloc[-1]),
                      int(s.n_variants.iloc[-1]))

    # reference line stops before the label column so it cannot strike a label
    ax.plot([35, 305], [0.75, 0.75], color=S.GREY, ls=":", lw=0.6, zorder=1)
    ax.text(38, 0.765, "75%", fontsize=5.5, color=S.GREY, va="bottom", ha="left")
    ax.set_xlabel("Cells per half")
    ax.set_ylabel("Fraction detectable")
    ax.set_xticks([50, 100, 150, 300])
    # Every plotted value lies between 0.507 and 1.000, so a 0-1.08 axis left the
    # lower half of the panel empty. Line plots do not require a zero baseline;
    # the axis is clipped to the occupied range (plus margin) and the 75% guide
    # keeps the "plateau below 75%" claim readable.
    ax.set_ylim(0.45, 1.05)
    ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_xlim(35, 395)
    S.despine(ax)
    ax.tick_params(length=2.2)

    # ---- direct end-of-curve labels (replace the repeated gene legend) ----
    for gene, (x, y, n_end) in ends.items():
        col = S.GENE_COLORS[gene]
        if x >= 300:                      # curve reaches the right edge
            ax.text(x * 1.06, y + 0.020, gene, color=col, fontsize=6,
                    fontweight="bold", va="center", ha="left")
            ax.text(x * 1.06, y - 0.032, f"n = {n_end}", color=S.GREY,
                    fontsize=5.2, va="center", ha="left")
        else:                             # TP53 stops at 150 (no 300-cell rung)
            ax.text(x + 10, y + 0.040, gene, color=col, fontsize=6,
                    fontweight="bold", va="center", ha="left")
            ax.text(x + 10, y - 0.012, f"n = {n_end}", color=S.GREY,
                    fontsize=5.2, va="center", ha="left")

    S.save(fig, "fig3e_titration")
    print(pw.pivot_table(index="gene", columns="n_sub", values="frac_detectable")
          .reindex(S.GENE_ORDER).round(3).to_string())


if __name__ == "__main__":
    main()
