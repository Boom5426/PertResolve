"""Fig 2g: representative held-out TP53 variants (Ridge-esm).

Message: direction of the predicted allele effect is stable across variants
(pearson_delta ~0.74-0.82), yet the allele *ranking* score (PDS_cos) spans the
whole 0-1 range. Similar direction, very different allele ranking.

Data: D.exttheta(), method=='Ridge-esm', gene=='TP53', per variant.
A representative subset of 10 distinct variants (deduplicated by variant, one
row each) evenly spaced across the observed PDS_cos range is shown.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

N_SHOW = 10
GENE = "TP53"
METHOD = "Ridge-esm"


def select_variants():
    df = D.exttheta()
    sub = df[(df.method == METHOD) & (df.gene == GENE)].copy()
    # one representative row per distinct variant, sorted by ranking score
    sub = sub.sort_values("PDS_cos").drop_duplicates("variant", keep="first")
    idx = np.linspace(0, len(sub) - 1, N_SHOW).round().astype(int)
    return sub.iloc[idx].reset_index(drop=True)


def main():
    S.apply_rcparams()
    pick = select_variants()
    x = np.arange(len(pick))
    labels = pick["variant"].tolist()
    pdir = pick["pearson_delta"].to_numpy()
    pds = pick["PDS_cos"].to_numpy()

    col = S.GENE_COLORS[GENE]

    fig, ax = S.panel(64, 46)

    # two horizontal bands: top = direction (stable), bottom = ranking (spread)
    Y_DIR = 1.0
    Y_PDS = 0.0

    # --- TOP band: pearson_delta lollipops (rescaled to sit as a thin strip) ---
    # map pearson_delta -> small vertical offset around Y_DIR so the strip reads
    # as "all high, all similar". Use its own mini-axis on the right.
    dir_lo, dir_hi = 0.70, 0.85
    dir_h = 0.34  # visual height of the top strip
    y_dir = Y_DIR + (pdir - dir_lo) / (dir_hi - dir_lo) * dir_h

    ax.vlines(x, Y_DIR, y_dir, color=S.LIGHT_GREY, lw=0.7, zorder=1)
    ax.scatter(x, y_dir, s=11, color=col, edgecolors="white",
               linewidths=0.35, zorder=3, clip_on=False)

    # --- BOTTOM band: PDS_cos lollipops (0..1 spread) ---
    pds_h = 0.80
    y_pds = Y_PDS + pds * pds_h
    ax.vlines(x, Y_PDS, y_pds, color=S.LIGHT_GREY, lw=0.7, zorder=1)
    ax.scatter(x, y_pds, s=11, color=col, edgecolors="white",
               linewidths=0.35, zorder=3, clip_on=False)

    # chance line on the PDS row (0.5)
    y_chance = Y_PDS + 0.5 * pds_h
    ax.axhline(y_chance, xmin=0.02, xmax=0.98, color=S.GREY, lw=0.6,
               ls=(0, (3, 2)), zorder=2)
    ax.text(len(pick) - 0.5, y_chance, "chance", ha="right", va="bottom",
            fontsize=5, color=S.GREY)

    # ---- axis cosmetics -------------------------------------------------------
    ax.set_xlim(-0.6, len(pick) - 0.4)
    ax.set_ylim(Y_PDS - 0.04, Y_DIR + dir_h + 0.12)
    # keep the two band labels clear of the plot area
    ax.margins(x=0)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90, fontsize=5)
    ax.tick_params(axis="x", length=0, pad=1.5)

    # custom y ticks: PDS band 0/0.5/1, dir band 0.70/0.85
    yt = [Y_PDS, y_chance, Y_PDS + pds_h,
          Y_DIR, Y_DIR + dir_h]
    ytl = ["0", "0.5", "1", "0.70", "0.85"]
    ax.set_yticks(yt)
    ax.set_yticklabels(ytl, fontsize=5)
    ax.tick_params(axis="y", length=1.6, pad=1.5)

    S.despine(ax, keep=("left", "bottom"))

    # band labels on the left in axes-fraction coords, clear of the tick numbers
    ymin, ymax = Y_PDS - 0.04, Y_DIR + dir_h + 0.12
    def yfrac(yval):
        return (yval - ymin) / (ymax - ymin)
    ax.text(-0.135, yfrac(Y_PDS + pds_h / 2), "PDS$_{cos}$ (ranking)",
            transform=ax.transAxes, ha="center", va="center", fontsize=5.5,
            rotation=90, color=S.INK)
    ax.text(-0.135, yfrac(Y_DIR + dir_h / 2), "pearson$_\\Delta$ (direction)",
            transform=ax.transAxes, ha="center", va="center", fontsize=5.5,
            rotation=90, color=S.INK)

    # message annotation
    ax.set_title("Held-out TP53 variants: similar direction,\nvery different allele ranking",
                 fontsize=5.8, pad=3, loc="left", linespacing=1.1)

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig2g_pervariant"))

    # echo numbers plotted
    for _, r in pick.iterrows():
        print(f"{r.variant:>7s}  split={r.split}  "
              f"pearson_delta={r.pearson_delta:.3f}  PDS_cos={r.PDS_cos:.3f}")


if __name__ == "__main__":
    main()
