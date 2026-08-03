"""Fig. 5c — phase map using ceiling CI and the full 18-head model range."""
from __future__ import annotations

W_MM, H_MM = 48.0, 31.0

import fig5_common as S


def main() -> None:
    S.apply_style()
    oracle = S.read_csv("oracle_ceiling.csv").query("scope != 'ALL'").set_index("scope")
    heads = S.per_gene_head_means()

    fig, ax = S.panel(W_MM, H_MM)
    fig.subplots_adjust(left=0.235, right=0.965, bottom=0.29, top=0.95)

    lo, hi = 0.44, 0.84
    ax.plot([lo, hi], [lo, hi], color=S.MID_GREY, lw=0.7,
            ls=(0, (3.2, 2.2)), zorder=0)
    ax.axvline(0.5, color=S.LIGHT_GREY, lw=0.6, zorder=0)
    ax.axhline(0.5, color=S.LIGHT_GREY, lw=0.6, zorder=0)

    # Every gene sits below the diagonal, so the upper-left half of the panel is
    # empty. Labels go there with short leaders rather than beside their markers:
    # TP53 and KRAS have nearly the same ceiling, so adjacent labels collided and
    # TP53's ran off the left edge of a 48 mm panel.
    label_anchors = {
        "TP53": (0.452, 0.662, "left"),
        "KRAS": (0.452, 0.752, "left"),
        "GATA1": (0.575, 0.600, "left"),
        "JAK1": (0.690, 0.690, "left"),
    }
    for gene in S.GENE_ORDER:
        row = oracle.loc[gene]
        values = heads.loc[heads["gene"] == gene, "PDS_cos"].to_numpy()
        x = float(row["PDS_oracle"])
        y = float(values.max())
        color = S.GENE_COLORS[gene]
        ax.plot([x, x], [float(values.min()), float(values.max())],
                color=color, lw=1.0, alpha=0.55, zorder=2)
        ax.errorbar(
            [x], [y],
            xerr=[[x - float(row["ci_lo"])], [float(row["ci_hi"]) - x]],
            fmt="D", markersize=4.5, color=color, markerfacecolor="white",
            markeredgecolor=color, markeredgewidth=0.9, ecolor=color,
            elinewidth=0.85, capsize=1.4, zorder=5,
        )
        lx, ly, ha = label_anchors[gene]
        ax.annotate(
            gene, xy=(x, y), xytext=(lx, ly), fontsize=5.8, color=S.INK,
            ha=ha, va="center", zorder=6,
            arrowprops=dict(arrowstyle="-", lw=0.45, color=S.MID_GREY,
                            shrinkA=1.0, shrinkB=3.0),
        )


    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([0.5, 0.6, 0.7, 0.8])
    ax.set_yticks([0.5, 0.6, 0.7, 0.8])
    # "native depth" is stated in the caption; at 48 mm the full label was 43 mm
    # of text on a 36 mm axis and ran off the panel.
    ax.set_xlabel("Replicate ceiling (PDS)")
    ax.set_ylabel("Best model PDS")
    S.despine(ax)

    S.save(fig, "fig5c_phase_map")


if __name__ == "__main__":
    main()
