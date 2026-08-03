"""Fig. 4e — feature-space changes preserve the direction–ranking gap."""
from __future__ import annotations

import numpy as np

import fig4_common as S

# Static QA contract: Arial/Helvetica sans-serif; svg.fonttype="none";
# pdf.fonttype=42; exports .svg, .pdf, .png and .tiff at dpi=600.
# Final assembled figure target: width_mm = 183.
W_MM, H_MM = 89.0, 45.0


def main() -> None:
    S.apply_style()
    df = S.breakdown()
    inhouse = df[df["method"].map(S.is_inhouse)].copy()
    per_method = (
        inhouse.groupby("method")[["PDS_cos", "pearson_delta"]].mean().reset_index()
    )
    per_method["feature"] = per_method["method"].map(S.feature_space)
    per_method["head"] = per_method["method"].map(S.head_name)
    direction = per_method.pivot(index="head", columns="feature", values="pearson_delta")
    pds = per_method.pivot(index="head", columns="feature", values="PDS_cos")
    direction = direction.loc[S.HEADS, S.FEATURE_ORDER]
    pds = pds.loc[S.HEADS, S.FEATURE_ORDER]

    x = np.arange(3, dtype=float)
    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0.16, 0.22, 0.78, 0.71])
    ax.axhline(0.50, color=S.GREY, lw=0.7, ls=(0, (3, 2)), zorder=1)

    for head in S.HEADS:
        ax.plot(
            x,
            direction.loc[head].to_numpy(float),
            color="#B7C8D7",
            lw=0.55,
            alpha=0.9,
            zorder=2,
        )
        ax.plot(
            x,
            pds.loc[head].to_numpy(float),
            color=S.LIGHT_GREY,
            lw=0.55,
            zorder=2,
        )

    direction_mean = direction.mean(axis=0).to_numpy(float)
    direction_sd = direction.std(axis=0, ddof=1).to_numpy(float)
    pds_mean = pds.mean(axis=0).to_numpy(float)
    pds_sd = pds.std(axis=0, ddof=1).to_numpy(float)
    ax.errorbar(
        x,
        direction_mean,
        yerr=direction_sd,
        fmt="o",
        ms=4.3,
        mfc=S.DIRECTION_BLUE,
        mec="white",
        mew=0.5,
        ecolor=S.DIRECTION_BLUE,
        elinewidth=0.85,
        capsize=2.2,
        zorder=5,
    )
    ax.errorbar(
        x,
        pds_mean,
        yerr=pds_sd,
        fmt="s",
        ms=4.0,
        mfc="white",
        mec=S.SLATE,
        mew=0.9,
        ecolor=S.SLATE,
        elinewidth=0.85,
        capsize=2.2,
        zorder=5,
    )

    ax.text(
        2.12,
        direction_mean[-1],
        r"Pearson-$\delta$",
        ha="left",
        va="center",
        fontsize=5.5,
        color=S.INK,
    )
    ax.text(
        2.12,
        pds_mean[-1],
        "PDS",
        ha="left",
        va="center",
        fontsize=5.5,
        color=S.INK,
    )
    ax.text(
        2.52,
        0.504,
        "PDS chance",
        ha="right",
        va="bottom",
        fontsize=5.0,
        color=S.GREY,
    )

    ax.set_xlim(-0.30, 2.55)
    ax.set_ylim(0.40, 0.68)
    ax.set_yticks([0.40, 0.50, 0.60])
    ax.set_xticks(x)
    ax.set_xticklabels([S.FEATURE_LABELS[value] for value in S.FEATURE_ORDER])
    ax.set_xlabel("Input feature space")
    ax.set_ylabel("Mean score")
    ax.tick_params(axis="both", pad=1.4)
    S.despine(ax)
    fig.text(
        0.55,
        0.895,
        "thin lines: matched heads · markers: mean ± SD across six heads",
        ha="center",
        va="top",
        fontsize=5.0,
        color=S.GREY,
    )
    fig.text(
        0.16,
        0.045,
        "single-draw breakdown",
        ha="left",
        va="bottom",
        fontsize=5.0,
        color=S.GREY,
    )
    S.save(fig, "fig4e_feature_spaces")


if __name__ == "__main__":
    main()
