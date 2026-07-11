"""Figure 3d - Effect size and sampling depth jointly determine the window.

Data-direct. Source: native-depth per-perturbation rows (fig3_data), plotting
variant-to-WT effect size (x, log) against the split-half ratio D_self/D_null (y),
with marker area proportional to effective cells per variant and colour by gene.
Rebuts the "TP53/KRAS are just low-depth" reading: they carry many cells yet sit
near R=1, whereas JAK1 has fewer cells but a wide (low-R) window. No regression is
fit; the panel shows the joint landscape only.

Run:  python fig3d_landscape.py  ->  fig3d_landscape.pdf (+ .png)
"""
from __future__ import annotations

import numpy as np
from matplotlib.lines import Line2D

import fig3_data as D
import nm_style as S


def main() -> None:
    S.apply_rcparams()
    fig, ax = S.panel(62, 50)

    for gene in S.GENE_ORDER:
        n = D.native_rankability(gene)
        es = n["effect_size"].to_numpy(float)
        ratio = n["ratio"].to_numpy(float)
        ncell = n["n_cells"].to_numpy(float)
        good = np.isfinite(es) & np.isfinite(ratio) & (es > 0)
        sizes = 3 + 22 * (ncell[good] - ncell[good].min()) / (np.ptp(ncell[good]) + 1e-9)
        ax.scatter(es[good], ratio[good], s=sizes, color=S.GENE_COLORS[gene],
                   alpha=0.55, linewidths=0, zorder=3)

    ax.axhline(1.0, color=S.GREY, ls="--", lw=0.7, zorder=1)
    ax.set_xscale("log")
    ax.set_xlabel(r"Variant-to-WT effect size (energy $D_\mathrm{null}$)")
    ax.set_ylabel(r"$D_\mathrm{self}\,/\,D_\mathrm{null}$")
    ax.set_ylim(-0.03, 1.4)
    S.despine(ax)
    ax.tick_params(length=2.2)

    gene_handles = [Line2D([0], [0], marker="o", color="none", markeredgecolor="none",
                           markerfacecolor=S.GENE_COLORS[g], markersize=4, label=g)
                    for g in S.GENE_ORDER]
    leg1 = ax.legend(handles=gene_handles, loc="lower left", title="gene",
                     handletextpad=0.2, labelspacing=0.25, borderpad=0.2,
                     title_fontsize=5.5)
    leg1.get_title().set_color(S.GREY)
    ax.add_artist(leg1)
    size_handles = [Line2D([0], [0], marker="o", color="none", markeredgecolor="none",
                           markerfacecolor=S.GREY, markersize=ms, label=lab)
                    for ms, lab in [(2.2, "few"), (5.0, "many")]]
    ax.legend(handles=size_handles, loc="lower right", title="cells/variant",
              handletextpad=0.2, labelspacing=0.4, borderpad=0.2, title_fontsize=5.5,
              labelcolor=S.GREY)

    S.save(fig, "fig3d_landscape")
    print("drawn: effect-size vs ratio, size=cells; genes", S.GENE_ORDER)


if __name__ == "__main__":
    main()
