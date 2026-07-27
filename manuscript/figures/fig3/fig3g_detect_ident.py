"""Figure 3d (composite letter) - Detection vs identification among sibling alleles.

Data-direct. Source: results/canonical/pairwise_resolvability.csv (D.pairwise()).
Two paired horizontal bars per gene, expressed as %:
  - "sibling pairs resolvable"  = frac_pairs_resolvable
  - "variants identifiable"     = frac_identifiable

Committed values:
  gene    frac_pairs_resolvable   frac_identifiable   median_nn_dist  median_Dself
  TP53    0.000                   0.000               0.230           0.965
  KRAS    0.000                   0.000               0.375           0.977
  GATA1   0.038                   0.000               0.461           0.981
  JAK1    0.785                   0.154               0.197           0.275

Message: allele-allele identification is ~0 for TP53/KRAS/GATA1 and only high
for JAK1 (78.5% of sibling pairs resolvable, 15.4% of variants identifiable),
because JAK1 sibling separation (nn_dist) is small relative to its self-noise
floor (D_self ~0.28) whereas the other three genes sit near D_self ~0.97.

Layout notes (Nature Methods pass): drawn at its final placed size so on-page type
is ~6 pt (composite scale ~1.0); every bar sits on a light 0-100% track so an
honest 0% reads as a measured zero on the full range rather than as empty canvas;
the two bar styles are named once in a quiet header key, and gene identity is
carried by the coloured y-axis labels (no gene legend is repeated here).

Encoding note: both quantities here are sibling-versus-sibling identification
measures, so this panel deliberately does NOT use an outline/solid (open/filled)
code. Open versus filled is reserved for panel e, where it means detection (vs
wild type) versus identification (sibling vs sibling); reusing it here for a
looser/stricter contrast would make open/filled mean two things in one row.
Instead the two series are coded by tint of the same gene hue: pale = the looser
pair-level measure, full = the stricter variant-level measure. The key therefore
shows filled swatches only, and cannot be misread as panel e's marker key.

Run:  python fig3g_detect_ident.py  ->  fig3g_detect_ident.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

import matplotlib.patches as mpatches
from matplotlib.colors import to_hex, to_rgb

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

W_MM, H_MM = 95.7, 31.8
XMAX = 100.0
TINT = 0.42          # fraction of white blended into the looser-measure bars


def tint(color: str, amount: float = TINT) -> str:
    """Blend ``color`` toward white as an opaque hex, so a pale bar drawn over the
    light-grey 0-100 % track stays a flat colour instead of compositing with it."""
    r, g, b = to_rgb(color)
    return to_hex((r + (1.0 - r) * amount,
                   g + (1.0 - g) * amount,
                   b + (1.0 - b) * amount))


def main() -> None:
    S.apply_rcparams()
    pw = D.pairwise().set_index("gene")

    fig, ax = S.panel(W_MM, H_MM)
    genes = S.GENE_ORDER

    group_h = 0.70          # total vertical span of a gene's pair of bars
    bar_h = group_h / 2.0
    centers = list(range(len(genes)))

    for c, gene in zip(centers, genes):
        color = S.GENE_COLORS[gene]
        resolv = float(pw.loc[gene, "frac_pairs_resolvable"]) * 100.0
        ident = float(pw.loc[gene, "frac_identifiable"]) * 100.0

        y_res = c - bar_h / 2.0
        y_id = c + bar_h / 2.0

        # full-range tracks: a 0% bar is a measured zero on a visible 0-100% scale
        for y in (y_res, y_id):
            ax.barh(y, XMAX, height=bar_h * 0.86, color=S.LIGHT_GREY, alpha=0.45,
                    edgecolor="none", zorder=1)

        # looser measure (sibling pairs resolvable): pale tint of the gene hue
        ax.barh(y_res, resolv, height=bar_h * 0.86, facecolor=tint(color),
                edgecolor="none", zorder=3)
        # stricter measure (variants identifiable): full-strength gene hue
        ax.barh(y_id, ident, height=bar_h * 0.86, color=color, alpha=0.95,
                edgecolor="none", zorder=3)

        # value labels just past the bar end (or just past the origin when 0)
        for y, val in ((y_res, resolv), (y_id, ident)):
            lab = f"{val:.1f}%" if val not in (0.0, 100.0) else "0%"
            ax.text(val + 1.6, y, lab, va="center", ha="left",
                    fontsize=5.5, color=S.INK, zorder=5)

    ax.set_yticks(centers)
    ax.set_yticklabels(genes)
    for t, g in zip(ax.get_yticklabels(), genes):
        t.set_color(S.GENE_COLORS[g])
        t.set_fontweight("bold")
    ax.invert_yaxis()
    ax.set_xlim(0, XMAX)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Sibling-allele identification (%)")
    S.despine(ax, keep=("bottom",))
    ax.tick_params(length=2.2)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(genes) - 0.42, -0.88)

    # one quiet style key on a single header line (neutral grey: it names the bar
    # TINT, not the gene colour, which the y-axis labels already carry). Both
    # swatches are filled, so this key cannot be read as panel e's open/filled code.
    pale = mpatches.Patch(facecolor=tint(S.GREY), edgecolor="none",
                          label="sibling pairs resolvable")
    solid = mpatches.Patch(facecolor=S.GREY, edgecolor="none",
                           label="variants identifiable")
    leg = ax.legend(handles=[pale, solid], loc="lower left",
                    bbox_to_anchor=(0.0, 1.005), ncol=2, handlelength=1.0,
                    handleheight=0.9, columnspacing=1.4, handletextpad=0.35,
                    borderaxespad=0.0, borderpad=0.0, fontsize=5.5)
    for t in leg.get_texts():
        t.set_color(S.GREY)
    leg.set_zorder(6)

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig3g_detect_ident"))

    # ---- report numbers ----
    for g in genes:
        print(f"{g:6} resolvable={pw.loc[g,'frac_pairs_resolvable']*100:5.1f}%  "
              f"identifiable={pw.loc[g,'frac_identifiable']*100:5.1f}%  "
              f"nn_dist={pw.loc[g,'median_nn_dist']:.3f}  "
              f"D_self={pw.loc[g,'median_Dself']:.3f}")


if __name__ == "__main__":
    main()
