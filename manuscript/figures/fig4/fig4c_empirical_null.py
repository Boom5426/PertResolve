"""Fig. 4c — empirical permutation calibration of PDS."""
from __future__ import annotations

import numpy as np

import fig4_common as S

# Static QA contract: Arial/Helvetica sans-serif; svg.fonttype="none";
# pdf.fonttype=42; exports .svg, .pdf, .png and .tiff at dpi=600.
# Final assembled figure target: width_mm = 183.
W_MM, H_MM = 105.0, 37.0


def main() -> None:
    S.apply_style()
    df = S.read_csv("permutation_null_pds.csv")
    df = df[df["metric"] == "PDS_cos"].copy()
    observed = df["obs_pds"].to_numpy(float)
    null_low = float(df["null_mean"].min())
    null_high = float(df["null_mean"].max())
    n = len(df)
    n_exact = int(np.isclose(observed, 0.50, atol=1e-6).sum())
    k_exceed = int((df["p_above_obs"] < 0.05).sum())
    exceed = 100.0 * k_exceed / n
    low, high = S.wilson_interval(k_exceed, n)
    low, high = 100.0 * low, 100.0 * high

    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0.10, 0.22, 0.58, 0.665])
    ax.axvspan(null_low, null_high, color=S.LIGHT_GREY, alpha=0.95, zorder=0)
    ax.hist(
        observed,
        bins=25,
        color=S.SLATE,
        edgecolor="white",
        linewidth=0.35,
        zorder=2,
    )
    ax.axvline(0.50, color=S.INK, lw=0.75, ls=(0, (3, 2)), zorder=3)
    ax.set_xlim(0.35, 0.60)
    ax.set_ylim(0, max(180, n_exact + 12))
    ax.set_xticks([0.35, 0.40, 0.45, 0.50, 0.55, 0.60])
    ax.set_xlabel("Observed PDS (cosine)")
    ax.set_ylabel("Combinations")
    ax.tick_params(axis="both", pad=1.3)
    S.despine(ax)
    ax.annotate(
        f"{n_exact}/{n} exactly 0.50\n({100.0 * n_exact / n:.0f}%; ties)",
        xy=(0.502, n_exact - 3),
        xytext=(0.530, n_exact * 0.72),
        ha="left",
        va="center",
        fontsize=5.2,
        color=S.INK,
        arrowprops=dict(arrowstyle="-", lw=0.5, color=S.GREY),
    )
    ax.text(
        0.02,
        0.96,
        f"null-mean range\n{null_low:.3f}–{null_high:.3f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.0,
        color=S.GREY,
        linespacing=1.15,
    )

    bx = fig.add_axes([0.75, 0.32, 0.21, 0.435])
    bx.axvline(5.0, color=S.GREY, lw=0.7, ls=(0, (3, 2)), zorder=1)
    bx.plot([low, high], [0.5, 0.5], color=S.SLATE, lw=1.2, zorder=2)
    bx.plot([low, low], [0.42, 0.58], color=S.SLATE, lw=0.8)
    bx.plot([high, high], [0.42, 0.58], color=S.SLATE, lw=0.8)
    bx.plot(
        exceed,
        0.5,
        marker="D",
        ms=4.4,
        mfc=S.PDS_BLUE,
        mec="white",
        mew=0.5,
        zorder=4,
    )
    bx.text(
        exceed,
        0.72,
        f"{exceed:.1f}%\n({k_exceed}/{n})",
        ha="center",
        va="bottom",
        fontsize=5.4,
        color=S.INK,
        linespacing=1.1,
    )
    # Right-aligned at the axis edge: left-aligned from 5.05 this label ran
    # 0.22 mm past the 89 mm canvas and was clipped in the composite.
    bx.text(8.0, 0.10, "expected 5%", ha="right", va="bottom", fontsize=5.0, color=S.GREY)
    bx.set_xlim(0, 8)
    bx.set_ylim(0, 1)
    bx.set_xticks([0, 4, 8])
    bx.set_yticks([])
    bx.tick_params(axis="x", pad=1.2)
    S.despine(bx, keep=("bottom",))
    fig.text(
        0.855,
        0.165,
        "Exceed own\n95th null (%)",
        ha="center",
        va="top",
        fontsize=5.0,
        color=S.INK,
        linespacing=1.05,
    )

    fig.text(
        0.50,
        0.905,
        "340 method × split × gene combinations · right: 95% Wilson interval",
        ha="center",
        va="top",
        fontsize=5.1,
        color=S.GREY,
    )
    S.save(fig, "fig4c_empirical_null")


if __name__ == "__main__":
    main()
