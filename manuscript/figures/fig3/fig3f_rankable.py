"""Figure 3h (composite letter) - Variants detection-rankable at native depth.

Data-direct. Source: results/unrankable_canonical.json (rankable = 100 - un-rankable),
cross-checked against fig3_data.rankable_fraction. Values (95% CI):
  TP53 0% (0-0), KRAS 0% (0-0), GATA1 2.4% (0.8-4.3), JAK1 90% (75-100).
This is the allele-level summary only; the gene-level atlases (Replogle/Norman/
Adamson) belong in Fig 6 / Extended Data, not here.

Layout notes (Nature Methods pass): drawn at its final placed size; each bar sits on
a light 0-100% track so an honest 0% reads as a measured zero on the full range;
value labels sit just past the bar (or past its CI); gene identity comes from the
coloured tick labels, so no gene legend is repeated.

Orientation note: this panel is drawn as VERTICAL columns while panel d is drawn as
paired horizontal bars. Both are four-gene, 0-100 % summaries on a light full-range
track, and side by side in the same figure they were being read as the same panel
twice. They measure different things (d: sibling-vs-sibling identification; h:
detection-level rankability at native depth), so the orientation and the single
CI-bearing series here mark them as distinct objects; each axis names its own
quantity in full.

Run:  python fig3f_rankable.py  ->  fig3f_rankable.pdf (+ .png)
"""
from __future__ import annotations

import fig3_data as D
import nm_style as S

# H_MM is set so that, at the composite's 53.95 mm placed width, this panel's
# placed height matches panels f/g and row 3 stays bottom-aligned (see
# fig3_assemble.tex). Re-measure the PDF page size if the axis furniture changes.
W_MM, H_MM = 58.5, 65.85
YMAX = 100.0


def main() -> None:
    S.apply_rcparams()
    uc = D.unrankable_canonical()

    fig, ax = S.panel(W_MM, H_MM)
    genes = S.GENE_ORDER
    xs = range(len(genes))
    bar_w = 0.60
    for x, gene in zip(xs, genes):
        rank = 100.0 - uc[gene]["ur"]
        lo = 100.0 - uc[gene]["ci_hi"]
        hi = 100.0 - uc[gene]["ci_lo"]
        color = S.GENE_COLORS[gene]
        # full-range track: a 0% column is a measured zero on a visible 0-100% scale
        ax.bar(x, YMAX, width=bar_w, color=S.LIGHT_GREY, alpha=0.45,
               edgecolor="none", zorder=1)
        ax.bar(x, rank, width=bar_w, color=color, alpha=0.95, edgecolor="none",
               zorder=3)
        ax.plot([x, x], [lo, hi], color=S.INK, lw=0.9, zorder=4, solid_capstyle="round")
        lab = f"{rank:.1f}%" if rank not in (0.0, 100.0) else f"{rank:.0f}%"
        ax.text(x, max(rank, hi) + 2.5, lab, ha="center", va="bottom", fontsize=6,
                color=S.INK, zorder=5)

    ax.set_xticks(list(xs))
    ax.set_xticklabels(genes)
    for t, g in zip(ax.get_xticklabels(), genes):
        t.set_color(S.GENE_COLORS[g]); t.set_fontweight("bold")
    ax.set_ylim(0, 118)
    ax.set_xlim(-0.65, len(genes) - 0.35)
    ax.set_ylabel("Variants rankable at native depth (%)")
    ax.set_yticks([0, 25, 50, 75, 100])
    S.despine(ax, keep=("left",))
    ax.tick_params(length=2.2)
    ax.tick_params(axis="x", length=0)

    S.save(fig, "fig3f_rankable")
    for g in genes:
        print(f"{g:6} rankable={100-uc[g]['ur']:.1f}%  n={uc[g]['n']}")


if __name__ == "__main__":
    main()
