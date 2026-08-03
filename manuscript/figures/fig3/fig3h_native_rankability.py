"""Fig. 3h — detection-rankable variants at native depth with Wilson CIs."""
from __future__ import annotations

import fig3_common as S

# Publication contract inherited from fig3_common: Arial/Helvetica sans-serif;
# svg.fonttype="none"; pdf.fonttype=42; exports .svg, .pdf, .png and .tiff at dpi=600.
# Final assembled figure target: width_mm = 183.
W_MM, H_MM = 87.0, 36.0


def main() -> None:
    S.apply_style()
    native = S.read_json("unrankable_canonical.json")["canonical_native"]
    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0.23, 0.24, 0.73, 0.57])

    counts = {}
    for gene in S.GENE_ORDER:
        n = int(native[gene]["n"])
        rank_pct = 100.0 - float(native[gene]["ur"])
        k = int(round(rank_pct * n / 100.0))
        counts[gene] = (k, n)

    y_positions = list(range(len(S.GENE_ORDER)))[::-1]
    for y, gene in zip(y_positions, S.GENE_ORDER):
        k, n = counts[gene]
        estimate = 100.0 * k / n
        lo, hi = S.wilson_interval(k, n)
        lo, hi = 100.0 * lo, 100.0 * hi
        color = S.GENE_COLORS[gene]
        ax.plot([0, 100], [y, y], color=S.PALE_GREY, lw=0.55, zorder=0)
        ax.plot([lo, hi], [y, y], color=color, lw=1.1, zorder=2)
        ax.plot([lo, lo], [y - 0.10, y + 0.10], color=color, lw=0.8)
        ax.plot([hi, hi], [y - 0.10, y + 0.10], color=color, lw=0.8)
        ax.plot(
            estimate,
            y,
            marker="o",
            ms=4.6,
            mfc=color,
            mec="white",
            mew=0.6,
            zorder=4,
            clip_on=False,
        )
        label_x = max(estimate, hi) + 2.0
        label = f"{estimate:.1f}%" if estimate not in (0.0, 100.0) else f"{estimate:.0f}%"
        ax.text(
            label_x,
            y,
            label,
            ha="left",
            va="center",
            fontsize=6.0,
            color=S.INK,
        )

    labels = [f"{g}\n{counts[g][0]}/{counts[g][1]}" for g in S.GENE_ORDER]
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, linespacing=1.08)
    S.gene_ticklabels(ax, "y")
    ax.set_xlim(-4, 106)
    ax.set_ylim(-0.45, 3.45)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Detection-rankable variants (%)", labelpad=2)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", pad=1.5)
    S.despine(ax, keep=("bottom",))
    fig.text(
        0.23,
        0.855,
        "point: proportion · line: 95% Wilson CI",
        ha="left",
        va="center",
        fontsize=5.0,
        color=S.GREY,
    )
    fig.text(
        0.595,
        0.965,
        r"Native depth · primary criterion $S>W$",
        ha="center",
        va="top",
        fontsize=6.1,
    )
    S.save(fig, "fig3h_native_rankability")


if __name__ == "__main__":
    main()
