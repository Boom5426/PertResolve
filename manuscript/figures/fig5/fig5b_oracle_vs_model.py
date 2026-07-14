"""Figure 5b - Per-gene oracle ceiling vs best honest model, against the 0.5 chance line.

Data-direct. Two committed sources:
  - Oracle ceiling (filled marker + bootstrap-CI whisker): results/_remote/unified/
    oracle_ceiling.csv, native-depth per-gene PDS_oracle (TP53 0.485 [0.445,0.525],
    KRAS 0.500 [0.440,0.566], GATA1 0.572 [0.549,0.595], JAK1 0.792 [0.742,0.838]).
  - Best honest model (open marker): max over the 18 in-house feature-model heads of
    the per-gene mean PDS_cos in results/canonical/unified_results5.csv (canonical
    multi-seed harness; TP53 0.51, KRAS 0.54, GATA1 0.50, JAK1 0.52).

Message: TP53/KRAS oracle ~= chance (measurement-limited, no headroom to compete for);
GATA1 intermediate; JAK1 oracle 0.79 >> best model 0.48 (real computational gap). The
best honest model sits near chance for every gene; only JAK1 leaves room to improve.

Run:  python fig5b_oracle_vs_model.py  ->  fig5b_oracle_vs_model.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

CHANCE = 0.50
REGIME = {
    "TP53": "measurement-\nlimited",
    "KRAS": "measurement-\nlimited",
    "GATA1": "intermediate",
    "JAK1": "computational\ngap",
}


def best_model_pds() -> dict[str, float]:
    """Per-gene best in-house head, canonical multi-seed harness (D.best_model)."""
    return D.best_model()


def main() -> None:
    S.apply_rcparams()
    orc = D.oracle().set_index("scope")
    best = best_model_pds()

    fig, ax = S.panel(58, 46)

    genes = S.GENE_ORDER  # TP53, KRAS, GATA1, JAK1
    ypos = {g: i for i, g in enumerate(genes)}  # 0..3 bottom->top

    # chance reference
    ax.axvline(CHANCE, color=S.GREY, lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax.text(CHANCE, len(genes) - 0.30, "chance", fontsize=5, color=S.GREY,
            ha="center", va="bottom")

    printed = []
    for g in genes:
        y = ypos[g]
        col = S.GENE_COLORS[g]
        o = float(orc.loc[g, "PDS_oracle"])
        lo = float(orc.loc[g, "ci_lo"])
        hi = float(orc.loc[g, "ci_hi"])
        b = float(best[g])

        # connector between best model and oracle
        ax.plot([min(b, lo), max(b, hi)], [y, y], "-", color=S.LIGHT_GREY,
                lw=1.0, zorder=2, solid_capstyle="round")
        # oracle CI whisker
        ax.plot([lo, hi], [y, y], "-", color=col, lw=1.0, zorder=3,
                solid_capstyle="round")
        # oracle: filled marker
        ax.scatter([o], [y], s=26, marker="o", color=col, edgecolor="white",
                   linewidths=0.5, zorder=5)
        # best model: open marker
        ax.scatter([b], [y], s=22, marker="o", facecolor="white", edgecolor=col,
                   linewidths=0.9, zorder=5)
        printed.append((g, b, o, lo, hi))

    # right-margin regime labels
    for g in genes:
        y = ypos[g]
        ax.text(1.005, y, REGIME[g], transform=ax.get_yaxis_transform(),
                fontsize=5, color=S.GENE_COLORS[g], ha="left", va="center",
                linespacing=0.95)

    # in-panel legend (filled=oracle, open=model), drawn as proxy markers
    lx = 0.585
    ax.scatter([lx], [0.16], s=26, marker="o", color=S.INK, edgecolor="white",
               linewidths=0.5, zorder=6)
    ax.text(lx + 0.012, 0.16, "oracle ceiling", fontsize=5, color=S.INK,
            ha="left", va="center")
    ax.scatter([lx], [-0.16], s=22, marker="o", facecolor="white",
               edgecolor=S.INK, linewidths=0.9, zorder=6)
    ax.text(lx + 0.012, -0.16, "best in-house model", fontsize=5, color=S.INK,
            ha="left", va="center")

    ax.set_yticks(list(ypos.values()))
    ax.set_yticklabels(genes)
    for lab in ax.get_yticklabels():
        lab.set_color(S.GENE_COLORS[lab.get_text()])
    ax.set_ylim(-0.55, len(genes) - 0.30)
    ax.set_xlim(0.42, 0.86)
    ax.set_xticks([0.5, 0.6, 0.7, 0.8])
    ax.set_xlabel("PDS (perturbation discrimination score)")
    S.despine(ax)
    ax.tick_params(length=2.2)

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig5b_oracle_vs_model"))

    print("gene  best_model  oracle  ci_lo  ci_hi")
    for g, b, o, lo, hi in printed:
        print(f"{g:6s} {b:.3f}      {o:.3f}  {lo:.3f}  {hi:.3f}")


if __name__ == "__main__":
    main()
