"""Fig 6f: native vs depth-matched (n=50/half) un-rankable fraction on external
gene-level atlases (Replogle, Adamson, Norman).

Message: depth normalization raises Adamson and Norman toward Replogle, and
Replogle stays high -> a dataset-specific effect-size structure, not depth alone.

Pairing correction (2026-07-31). The panel used to read both arms from
``D.canonical()``, which computes each arm over every perturbation that HAS that
arm. Those sets are not the same: at native depth Replogle has 1,832
perturbations and Norman 236, while only 1,299 and 224 respectively reach a
50-cell-per-half rung. Drawing a connector between two differently-denominated
fractions asserts a within-perturbation shift that was never computed, and for
Replogle it inverts the sign: unpaired the fraction appears to fall 55.3 -> 54.9,
while on the 1,299 perturbations present in both arms it RISES 52.0 -> 54.9.

Both arms are therefore computed here on the paired subset, from the committed
``results/second_probe_rankability_table.csv``. The unpaired path is recomputed
too and asserted against ``results/canonical_numbers.json``, so this panel cannot
drift away from the canonical file without failing loudly.

Adamson is unaffected (96 perturbations in both arms).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import fig6_data as D   # in the same fig6/ dir

import numpy as np
import pandas as pd

MATCH_N = 50          # cells per split half in the depth-matched arm
METRIC, SPACE = "edist", "pca"
DATASETS = ["Replogle", "Adamson", "Norman"]   # top-to-bottom display order


def _table() -> pd.DataFrame:
    t = pd.read_csv(S.repo_root() / "results" / "second_probe_rankability_table.csv")
    return t[(t.metric == METRIC) & (t.space == SPACE) & (t.n_cells > 0)]


def arms(ds: str, tab: pd.DataFrame) -> dict:
    """Un-rankable percentages for one dataset, paired and unpaired.

    native  = each perturbation at its deepest available split-half rung
    matched = the same perturbation at the 50-cell-per-half rung
    """
    sub = tab[tab.dataset == ds]
    native = sub.loc[sub.groupby("perturbation")["n_work"].idxmax()]
    matched = sub[sub.n_work == MATCH_N].groupby("perturbation")["rankable"].max()
    both = native.perturbation.isin(matched.index)

    def pct(rankable: np.ndarray) -> tuple[float, int, int]:
        n = len(rankable)
        k = int(n - np.sum(rankable))        # un-rankable count
        return 100.0 * k / n, k, n

    p_nat, k_nat, n_nat = pct(native.loc[both, "rankable"].to_numpy())
    p_mat, k_mat, n_mat = pct(matched.reindex(native.loc[both, "perturbation"]).to_numpy())
    assert n_nat == n_mat, f"{ds}: paired arms disagree on n ({n_nat} vs {n_mat})"
    return dict(
        native=p_nat, native_k=k_nat, matched=p_mat, matched_k=k_mat, n=n_nat,
        # unpaired, for the canonical-file guard only; never plotted
        native_unpaired=pct(native["rankable"].to_numpy())[0],
        matched_unpaired=pct(matched.to_numpy())[0],
        n_native_unpaired=len(native), n_matched_unpaired=len(matched),
    )


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion, in percent."""
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    denom = 1.0 + z * z / n
    centre = p + z * z / (2 * n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return 100.0 * (centre - half) / denom, 100.0 * (centre + half) / denom


def main() -> None:
    S.apply_rcparams()
    tab = _table()
    A = {d: arms(d, tab) for d in DATASETS}

    # Guard: the unpaired recomputation must still equal the canonical file.
    c = D.canonical()
    for d in DATASETS:
        for key, ref in (("native_unpaired", c["unrankable_native_pct"][d]),
                         ("matched_unpaired", c["unrankable_matched50_pct"][d])):
            got = A[d][key]
            assert abs(got - ref) < 0.05, (
                f"{d} {key}: recomputed {got:.2f}% but canonical_numbers.json says "
                f"{ref:.2f}%. The table and the canonical file have diverged; fix "
                f"the source of truth before redrawing this panel.")

    y_pos = {d: i for i, d in enumerate(reversed(DATASETS))}
    fig, ax = S.panel(52.0, 44.0)

    for d in DATASETS:
        y = y_pos[d]
        col = D.DATASET_COLORS[d]
        xn, xm = A[d]["native"], A[d]["matched"]
        nlo, nhi = wilson(A[d]["native_k"], A[d]["n"])
        mlo, mhi = wilson(A[d]["matched_k"], A[d]["n"])

        # Wilson 95% intervals first, so the markers sit on top of them
        ax.plot([nlo, nhi], [y - 0.16, y - 0.16], color=col, lw=0.7, alpha=0.55,
                solid_capstyle="butt", zorder=1)
        ax.plot([mlo, mhi], [y + 0.16, y + 0.16], color=col, lw=0.7, alpha=0.55,
                solid_capstyle="butt", zorder=1)
        # connecting dumbbell: legitimate now that both arms share one denominator
        ax.plot([xn, xm], [y, y], color=col, lw=0.9, zorder=2,
                solid_capstyle="round")
        ax.plot(xm, y, marker="o", ms=4.6, mfc="white", mec=col, mew=1.0,
                linestyle="none", zorder=3)
        ax.plot(xn, y, marker="o", ms=4.0, mfc=col, mec="white", mew=0.5,
                linestyle="none", zorder=4)

        # value labels: lower value left of its marker, higher value right
        lo, hi = min(xn, xm), max(xn, xm)
        if lo < 8.0:   # near the axis: lift above so it clears the tick label
            ax.text(lo, y + 0.24, f"{lo:.1f}", ha="center", va="bottom",
                    fontsize=5.5, color=S.INK)
        else:
            ax.text(lo - 2.6, y, f"{lo:.1f}", ha="right", va="center",
                    fontsize=5.5, color=S.INK)
        ax.text(hi + 2.6, y, f"{hi:.1f}", ha="left", va="center",
                fontsize=5.5, color=S.INK)

    ax.set_xlim(0, 60)
    ax.set_xticks([0, 20, 40, 60])
    ax.set_ylim(-0.5, len(DATASETS) - 0.5)
    ax.set_yticks(list(y_pos.values()))
    # n is part of the row label: the two arms are paired, so one n serves both
    ax.set_yticklabels([f"{d}\nn = {A[d]['n']:,}" for d in reversed(DATASETS)],
                       fontsize=6.5, linespacing=1.15)
    ax.set_xlabel("un-rankable perturbations (%)", fontsize=6.5)
    ax.tick_params(axis="both", length=2.0, pad=1.5)

    S.despine(ax, keep=("left", "bottom"))
    for tick, d in zip(ax.get_yticklabels(), reversed(DATASETS)):
        tick.set_color(D.DATASET_COLORS[d])

    # LEGEND ECONOMY: no detached legend. Direct-labelled once, above the Adamson
    # row, whose two markers are furthest apart.
    ax.text(A["Adamson"]["native"], 1.26, "native", ha="center", va="bottom",
            fontsize=5.5, color=S.INK)
    ax.text(A["Adamson"]["matched"], 1.26, "depth-matched\n(n = 50 per half)",
            ha="center", va="bottom", fontsize=5.5, color=S.INK, linespacing=1.2)

    fig.subplots_adjust(left=0.26, right=0.97, top=0.94, bottom=0.16)
    stem = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig6f_depthmatch")
    S.save(fig, stem)

    for d in DATASETS:
        a = A[d]
        print(f"{d:9s} n={a['n']:5d}  native {a['native']:5.1f}% "
              f"[{wilson(a['native_k'], a['n'])[0]:.1f},{wilson(a['native_k'], a['n'])[1]:.1f}]"
              f"  matched {a['matched']:5.1f}% "
              f"[{wilson(a['matched_k'], a['n'])[0]:.1f},{wilson(a['matched_k'], a['n'])[1]:.1f}]"
              f"   (unpaired was {a['native_unpaired']:.1f} -> {a['matched_unpaired']:.1f}"
              f", n {a['n_native_unpaired']} vs {a['n_matched_unpaired']})")


if __name__ == "__main__":
    main()
