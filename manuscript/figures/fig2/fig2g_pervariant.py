"""Fig. 2g — all distinct TP53 variants for Ridge ESM.

All 59 distinct variants are included. When a variant is evaluated in more than
one eligible split, scores are averaged across those rows before plotting.
Variants are ordered by residue position, not by either displayed outcome.
The identical 0–1 y scales and IQR bands show stable direction but widely
variable allele discrimination without outcome-driven example selection.

Run: python fig2g_pervariant.py

Static QA contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42; outputs
.svg, .pdf, .png and .tiff at dpi=600. Final assembled figure target:
width_mm = 183.
"""
from __future__ import annotations

import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

HERE = os.path.dirname(os.path.abspath(__file__))
GENE = "TP53"
METHOD = "Ridge-esm"

W_MM, H_MM = 107.5, 44.0
# Absolute millimetres, so they track the canvas width rather than a
# fraction of it: 102.5 was sized for the old 120.5 mm panel and left the
# 393 tick 11 mm outside the narrowed one.
AX_LEFT_MM, AX_WIDTH_MM = 15.0, 89.5
AX_H_MM = 12.5
AX_BOT_MM = 7.8
AX_TOP_MM = 25.2


def residue_position(name: str) -> int:
    match = re.search(r"(\d+)", str(name))
    if not match:
        raise ValueError(f"variant has no residue position: {name}")
    return int(match.group(1))


def aggregate_variants():
    df = D.exttheta()
    sub = df[(df.method == METHOD) & (df.gene == GENE)].copy()
    agg = (sub.groupby("variant", as_index=False)
              .agg(pearson_delta=("pearson_delta", "mean"),
                   PDS_cos=("PDS_cos", "mean"),
                   n_split_rows=("split", "size")))
    agg["position"] = agg.variant.map(residue_position)
    agg = agg.sort_values(["position", "variant"], kind="stable").reset_index(drop=True)

    # Small symmetric offsets reveal multiple substitutions at the same residue
    # without changing their residue-order interpretation.
    agg["x"] = agg["position"].astype(float)
    for _, idx in agg.groupby("position").groups.items():
        idx = list(idx)
        if len(idx) > 1:
            agg.loc[idx, "x"] += np.linspace(-1.5, 1.5, len(idx))
    return agg


def canvas():
    fig = plt.figure(figsize=(W_MM * S.MM, H_MM * S.MM))
    bg = fig.add_axes([0, 0, 1, 1])
    bg.set_axis_off()
    bg.set_xlim(0, 1)
    bg.set_ylim(0, 1)
    bg.add_patch(Rectangle((0, 0), 1, 1, facecolor="white",
                           edgecolor="none", zorder=-10))

    def band(bottom):
        return fig.add_axes([AX_LEFT_MM / W_MM, bottom / H_MM,
                             AX_WIDTH_MM / W_MM, AX_H_MM / H_MM])

    return fig, band(AX_TOP_MM), band(AX_BOT_MM)


def draw_distribution(ax, x, values, color, *, show_chance=False):
    q1, med, q3 = np.percentile(values, [25, 50, 75])
    ax.axhspan(q1, q3, color=color, alpha=0.10, zorder=0)
    ax.axhline(med, color=color, lw=0.65, alpha=0.75, zorder=1)
    if show_chance:
        ax.axhline(0.50, color=S.GREY, lw=0.6, ls=(0, (3, 2)), zorder=1)
        ax.text(389, 0.515, "chance", ha="right", va="bottom",
                fontsize=5.0, color=S.GREY)
    ax.scatter(x, values, s=9.5, facecolor=color, edgecolor="white",
               linewidth=0.3, alpha=0.82, zorder=3)
    label_y = min(0.97, med + 0.055)
    label_va = "bottom"
    if show_chance and med < 0.50:
        label_y = max(0.03, med - 0.055)
        label_va = "top"
    ax.text(389, label_y, f"median {med:.2f}",
            ha="right", va=label_va, fontsize=5.0, color=color,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.3), zorder=5)


def main() -> None:
    S.apply_rcparams()
    data = aggregate_variants()
    assert len(data) == 59

    color = S.GENE_COLORS[GENE]
    fig, ax_dir, ax_pds = canvas()
    draw_distribution(ax_dir, data.x, data.pearson_delta, color)
    draw_distribution(ax_pds, data.x, data.PDS_cos, color, show_chance=True)

    for ax in (ax_dir, ax_pds):
        ax.set_xlim(0, S.PROTEIN_LEN[GENE])
        ax.set_ylim(0, 1.03)
        ax.set_yticks([0, 0.5, 1.0])
        ax.tick_params(axis="y", length=1.8, pad=1.3, labelsize=5.5)
        S.despine(ax, keep=("left", "bottom"))
        ax.spines["left"].set_bounds(0, 1)

    ax_dir.set_xticks([])
    ax_dir.set_ylabel(r"Pearson-$\delta$", labelpad=1.4)
    ax_pds.set_xticks([0, 100, 200, 300, 393])
    ax_pds.tick_params(axis="x", length=1.8, pad=1.3, labelsize=5.5)
    ax_pds.set_ylabel(r"PDS$_{cos}$", labelpad=1.4)
    ax_pds.set_xlabel("TP53 residue position (Ridge ESM)", labelpad=1.4)

    fig.text(0.985, 0.965,
             "all 59 distinct variants · repeated split rows averaged per variant",
             ha="right", va="top", fontsize=5.0, color=S.GREY)

    S.save(fig, os.path.join(HERE, "fig2g_pervariant"), exact=True,
           formats=("pdf", "png", "svg", "tiff"))


if __name__ == "__main__":
    main()
