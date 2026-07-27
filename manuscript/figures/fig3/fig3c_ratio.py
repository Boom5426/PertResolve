"""Figure 3c (hero) - Per-variant split-half ratios define gene-specific regimes.

Data-direct. Source: results/second_probe_rankability_table.csv via the canonical
native-depth aggregation (fig3_data.native_rankability), verified to reproduce
canonical_numbers.json; the aggregate mean + 95% bootstrap CI overlay comes from
results/bootstrap_CIs.json. Per-gene D_self/D_null (mean, CI):
  TP53 0.96 (0.95-0.98), KRAS 1.00 (0.99-1.02), GATA1 0.88 (0.86-0.90), JAK1 0.21 (0.12-0.33).
Message: replicate noise approaches the variant signal for TP53/KRAS/GATA1 but not JAK1.

Run:  python fig3c_ratio.py   ->  fig3c_ratio.pdf (+ .png)
"""
from __future__ import annotations

import numpy as np

import fig3_data as D
import nm_style as S

RNG = np.random.default_rng(0)


def main() -> None:
    S.apply_rcparams()
    ci = D.dself_dnull_ci()

    fig, ax = S.panel(67.3, 64.4)

    # soft, annotation-only regime bands (NOT hard thresholds)
    ax.axhspan(0.90, 1.25, color=S.LIGHT_GREY, alpha=0.35, zorder=0)  # near floor
    ax.axhspan(0.0, 0.50, color=S.LIGHT_GREY, alpha=0.13, zorder=0)  # wide window
    ax.text(3.62, 1.14, "near floor", fontsize=5.5, color=S.GREY, ha="right", va="center")
    ax.text(-0.72, 0.30, "wide window", fontsize=5.5, color=S.GREY, ha="left", va="center")

    for i, gene in enumerate(S.GENE_ORDER):
        r = D.native_rankability(gene)["ratio"].to_numpy()
        r = r[np.isfinite(r)]
        color = S.GENE_COLORS[gene]

        parts = ax.violinplot(r, positions=[i], widths=0.85, showextrema=False)
        for b in parts["bodies"]:
            v = b.get_paths()[0].vertices
            v[:, 0] = np.clip(v[:, 0], i, np.inf)   # right half only
            b.set_facecolor(color); b.set_edgecolor("none"); b.set_alpha(0.30)

        jit = i - 0.06 - RNG.uniform(0, 0.24, size=r.size)
        ax.scatter(jit, r, s=2.2, color=color, alpha=0.5, linewidths=0, zorder=3)

        # canonical mean + 95% bootstrap CI (quoted values)
        m, lo, hi = ci[gene]["mean"], ci[gene]["lo"], ci[gene]["hi"]
        ax.plot([i + 0.02, i + 0.02], [lo, hi], color=S.INK, lw=1.0, zorder=5,
                solid_capstyle="round")
        ax.plot([i - 0.30, i + 0.30], [m, m], color=S.INK, lw=1.2, zorder=6)
        ax.text(i + 0.34, m, f"{m:.2f}", fontsize=6, va="center", color=S.INK)

    ax.axhline(1.0, color=S.GREY, ls="--", lw=0.7, zorder=1)
    ax.text(-0.72, 1.0, "R = 1", fontsize=5.5, color=S.GREY, va="bottom", ha="left")
    ax.set_xticks(range(4))
    ax.set_xticklabels(S.GENE_ORDER)
    for t, g in zip(ax.get_xticklabels(), S.GENE_ORDER):
        t.set_color(S.GENE_COLORS[g]); t.set_fontweight("bold")
    ax.set_xlim(-0.78, 3.7)
    ax.set_ylim(-0.03, 1.35)
    ax.set_ylabel(r"$D_\mathrm{self}/D_\mathrm{null}$")
    S.despine(ax)
    ax.tick_params(length=2.2)

    S.save(fig, "fig3c_ratio")
    for g in S.GENE_ORDER:
        n = D.native_rankability(g)
        print(f"{g:6} n={len(n):3d} ratio mean={n['ratio'].mean():.3f} "
              f"CI[{ci[g]['lo']:.3f},{ci[g]['hi']:.3f}]")


if __name__ == "__main__":
    main()
