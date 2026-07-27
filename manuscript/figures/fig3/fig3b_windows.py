"""Figure 3b - Representative closed and open measurement windows.

Data-direct. Source: results/_remote/unified/fig3b_selfnull_dist.csv, computed on the
server from raw cells (grid_cbv.npz) in the canonical PCA-50 space: 200 split-half seeds
per representative variant, giving the within-variant replicate distance D_self
(half 1 vs half 2, energy distance) and the variant-to-WT signal D_null.
Representative variants (nearest to each gene's median ratio):
  TP53 P222P  (closed window): D_self ~ 0.26 ~= D_null ~ 0.28  -> ratio ~ 0.96
  JAK1 R108Q  (open window):   D_self ~ 0.72 <<  D_null ~ 4.52 -> ratio ~ 0.16
Two sub-panels (own y-scale each) because the energy distances differ ~15x between genes.

Run:  python fig3b_windows.py  ->  fig3b_windows.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S


def _halfviolin(ax, x, data, color, side=1):
    parts = ax.violinplot(data, positions=[x], widths=0.8, showextrema=False)
    for b in parts["bodies"]:
        v = b.get_paths()[0].vertices
        v[:, 0] = np.clip(v[:, 0], x, np.inf) if side > 0 else np.clip(v[:, 0], -np.inf, x)
        b.set_facecolor(color); b.set_edgecolor("none"); b.set_alpha(0.55)
    med = np.median(data)
    ax.plot([x - 0.28, x + 0.28], [med, med], color=S.INK, lw=1.1, zorder=5)
    return med


def main() -> None:
    S.apply_rcparams()
    df = pd.read_csv(S.repo_root() / "results" / "canonical" / "fig3b_selfnull_dist.csv")

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(54.8 * S.MM, 55.0 * S.MM))

    panels = [
        ("TP53", axes[0], "closed window", S.GENE_COLORS["TP53"]),
        ("JAK1", axes[1], "open window", S.GENE_COLORS["JAK1"]),
    ]
    for gene, ax, regime, gcol in panels:
        d = df[df.gene == gene]
        var = d.variant.iloc[0]
        m_self = _halfviolin(ax, 0, d.D_self.to_numpy(), gcol)
        m_null = _halfviolin(ax, 1, d.D_null.to_numpy(), S.GREY)
        ratio = m_self / m_null
        ax.set_xticks([0, 1])
        ax.set_xticklabels([r"$D_\mathrm{self}$" + "\n(noise)", r"$D_\mathrm{null}$" + "\n(signal)"],
                           fontsize=6.0)
        ax.set_xlim(-0.6, 1.6)
        top = float(d.D_null.max()) * 1.15
        ax.set_ylim(0, top)
        ax.set_title(f"{gene} {var}\n{regime}  ($R$ = {ratio:.2f})", fontsize=6,
                     color=gcol, pad=3)
        S.despine(ax)
        ax.tick_params(length=2.2)
        # both sub-panels carry the y label: their scales differ ~15x, so an
        # unlabelled second axis invites a direct visual comparison that is wrong
        ax.set_ylabel("Energy distance")

    fig.subplots_adjust(wspace=0.5)
    S.save(fig, "fig3b_windows")
    print(df.groupby("gene")[["D_self", "D_null"]].median().round(3).to_string())


if __name__ == "__main__":
    main()
