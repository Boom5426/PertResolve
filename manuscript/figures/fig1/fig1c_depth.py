"""Figure 1c - Per-variant cell depth spans distinct sampling regimes.

Data-direct panel. Source: data/allele_perturb_bench.csv (committed), column
`n_cells` per variant. Reproduces the manuscript Fig 1c numbers exactly:
  median cells/variant  TP53 929, KRAS 1000, GATA1 355, JAK1 104
  total cells/gene      TP53 83.4k, KRAS 83.6k, GATA1 149.2k, JAK1 4.9k
  variants              98 / 93 / 255 / 26   (total 472; 321,043 cells)

Half-violin (density) + jittered points + median marker, log y. No result is
stated here; the depth spread only sets up the Fig 3 measurement analysis.

Run:  python fig1c_depth.py
Out:  fig1c_depth.pdf (+ .png preview)
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

RNG = np.random.default_rng(0)  # explicit seed for the point jitter only


def main() -> None:
    S.apply_rcparams()
    df = S.load_bench(exclude_wt=True)  # per-variant depth excludes the 2 WT rows

    fig, ax = S.panel(60, 41)
    for i, gene in enumerate(S.GENE_ORDER):
        vals = df.loc[df.gene == gene, "n_cells"].to_numpy(dtype=float)
        color = S.GENE_COLORS[gene]

        # half-violin (KDE) on the right side of each gene position
        parts = ax.violinplot(np.log10(vals), positions=[i], widths=0.9,
                              showextrema=False)
        for body in parts["bodies"]:
            verts = body.get_paths()[0].vertices
            verts[:, 0] = np.clip(verts[:, 0], i, np.inf)  # keep right half only
            body.set_facecolor(color)
            body.set_edgecolor("none")
            body.set_alpha(0.22)

        # jittered raw points on the left (kept vector: only ~470 points total)
        jit = i - 0.06 - RNG.uniform(0, 0.22, size=vals.size)
        ax.scatter(jit, np.log10(vals), s=2.0, color=color, alpha=0.55,
                   linewidths=0, zorder=3)

        # median marker + label
        med = np.median(vals)
        ax.plot([i - 0.34, i + 0.30], [np.log10(med)] * 2, color=S.INK, lw=1.0,
                zorder=4, solid_capstyle="round")
        ax.text(i + 0.02, np.log10(med) + 0.06, f"{med:.0f}", ha="center",
                va="bottom", fontsize=6, color=S.INK)
        # n variants below axis
        n = int((df.gene == gene).sum())
        ax.text(i, -0.14, f"n={n}", ha="center", va="top", fontsize=5.5,
                color=S.GREY, transform=ax.get_xaxis_transform())

    ax.set_xticks(range(len(S.GENE_ORDER)))
    # gene names in the gene colour: this panel is the figure's shared key for
    # "colour = gene" (direct labelling), so no detached gene legend is needed
    # anywhere in Figure 1.
    ax.set_xticklabels(S.GENE_ORDER)
    for lab, gene in zip(ax.get_xticklabels(), S.GENE_ORDER):
        lab.set_color(S.GENE_COLORS[gene])
        lab.set_fontweight("bold")
    ax.set_xlim(-0.6, len(S.GENE_ORDER) - 0.4)
    ax.set_ylabel("Cells per variant")

    # log-scale y with 10^k ticks shown as plain numbers; limits are tightened
    # to the observed range (33 to 1,987 cells per variant) so the panel is not
    # mostly empty axis
    ax.set_ylim(np.log10(25), np.log10(3000))
    ticks = [50, 100, 300, 1000, 3000]
    ax.set_yticks(np.log10(ticks))
    ax.set_yticklabels([f"{t:,}" for t in ticks])
    S.despine(ax)
    ax.tick_params(length=2.2)

    S.save(fig, "fig1c_depth")
    # report to stdout for the audit trail
    tot = df.groupby("gene")["n_cells"].agg(["median", "sum", "count"]).reindex(S.GENE_ORDER)
    print(tot.to_string())
    print("total cells:", int(df.n_cells.sum()), " total variants:", len(df))


if __name__ == "__main__":
    main()
