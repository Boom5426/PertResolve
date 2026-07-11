"""Figure 2f - Per-gene direction-vs-ranking dissociation gap.

Data-direct. Source: results/results_v4_exttheta.csv via remote_data.exttheta(),
restricted to the 18 in-house feature-model heads (remote_data.is_inhouse).
Per gene: gap = mean(pearson_delta) - mean(PDS_cos), pooled over the six splits
and all in-house predictors. Values:
  TP53 0.31, KRAS 0.19, GATA1 0.13, JAK1 -0.04.
Message: for TP53/KRAS/GATA1 predicted directions correlate well with truth yet
variant ranking (PDS_cos) lags far behind; for JAK1 the gap closes (~0), i.e. the
direction signal no longer outruns rankability.

Run:  python fig2f_gene_gap.py   ->  fig2f_gene_gap.pdf (+ .png)
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
    df = df[df["method"].apply(D.is_inhouse)]

    genes = S.GENE_ORDER  # TP53 > KRAS > GATA1 > JAK1
    gaps = {}
    for g in genes:
        sub = df[df["gene"] == g]
        gaps[g] = sub["pearson_delta"].mean() - sub["PDS_cos"].mean()

    fig, ax = S.panel(48, 44)

    ax.axhline(0.0, color=S.GREY, lw=0.6, zorder=1)

    xs = range(len(genes))
    for i, g in enumerate(genes):
        v = gaps[g]
        ax.bar(i, v, width=0.66, color=S.GENE_COLORS[g], edgecolor="none", zorder=3)
        # value label: above positive bars, below the (small) negative bar
        if v >= 0:
            ax.text(i, v + 0.012, f"{v:.2f}", ha="center", va="bottom",
                    fontsize=6, color=S.INK)
        else:
            ax.text(i, v - 0.012, f"{v:.2f}", ha="center", va="top",
                    fontsize=6, color=S.INK)

    ax.set_xticks(list(xs))
    ax.set_xticklabels(genes)
    for t, g in zip(ax.get_xticklabels(), genes):
        t.set_color(S.GENE_COLORS[g])
        t.set_fontweight("bold")

    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.set_ylim(-0.09, 0.37)
    ax.set_yticks([0.0, 0.1, 0.2, 0.3])
    ax.set_ylabel("direction minus ranking gap")

    # annotate the closing gap for JAK1 (message, not a threshold)
    ax.annotate("gap ~ 0", xy=(3, gaps["JAK1"]), xytext=(3, 0.13),
                ha="center", va="bottom", fontsize=5.5, color=S.GREY,
                arrowprops=dict(arrowstyle="-", lw=0.5, color=S.GREY,
                                shrinkA=0, shrinkB=2))

    S.despine(ax)
    ax.tick_params(length=2.2)

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig2f_gene_gap"))

    for g in genes:
        print(f"{g:6} gap={gaps[g]:+.4f}")


if __name__ == "__main__":
    main()
