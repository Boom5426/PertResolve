"""Fig 5c — measurement-limited vs computation-gap phase map.

Message: two distinct failure modes separate in ceiling-vs-model space.
  x = oracle discrimination ceiling per gene (independent second measurement).
  y = best honestly-scored model PDS per gene.

Source of truth (both axes committed):
  x: results/canonical/oracle_ceiling.csv via D.oracle() PDS_oracle
     (TP53 0.485, KRAS 0.500, GATA1 0.572, JAK1 0.792).
  y: results/canonical/unified_results5.csv via D.best_model() (best in-house head,
     canonical multi-seed harness): TP53 0.51, KRAS 0.54, GATA1 0.50, JAK1 0.52.
Both axes now regenerate from committed CSVs; the same y feeds Fig 5b and SI Table
tab:window "Best model PDS" so all three stay in sync.

Region shading is neutral grey only: blue/orange/purple/green are reserved for the
four genes throughout Figure 5, so no region may borrow a gene hue.

Regions (faint labels only, no leakage of ceiling/model values):
  diagonal    "benchmark-solvable"    ceiling above chance and model TRACKS the
              ceiling, i.e. a narrow band hugging y = x. The half-plane ABOVE y = x
              (model scores higher than its own oracle ceiling) is not an attainable
              regime: it is left unshaded and unlabelled, since any excursion into it
              is best-of-many-heads selection noise (see panel b).
  below band  "computation-limited"   ceiling clears chance, model does not follow.
  left of x = 0.5  "measurement-limited"  the ceiling itself is at or below chance.

Threshold correction (2026-07-31). The panel used to split the plane with a
vertical ``X_SEP = 0.55``, a constant that appeared nowhere in the caption, the
Results, the Methods or the SI. It decided a classification instead of
illustrating one. The only threshold drawn now is chance, 0.50, which is defined
throughout the paper, and every ceiling carries its committed 95% bootstrap
interval, so the reader can see directly which ceilings clear chance:

    TP53  0.485 [0.445, 0.525]   covers 0.50  -> ceiling not above chance
    KRAS  0.500 [0.440, 0.566]   covers 0.50  -> ceiling not above chance
    GATA1 0.572 [0.549, 0.595]   EXCLUDES 0.50
    JAK1  0.792 [0.742, 0.838]   EXCLUDES 0.50

Note what this shows and what it does not. On the oracle-PDS diagnostic GATA1's
ceiling is above chance, modestly, and its best model (0.50) does not follow it;
GATA1 therefore sits with JAK1 below the tracking band, not with TP53 and KRAS.
The Results call GATA1 measurement-limited on DIFFERENT diagnostics (D_self/D_null
0.88, 3.8% of sibling pairs resolvable, 2.4% of variants rankable at native
depth), which concern detection and identification rather than the attainable PDS.
The two statements are about different quantities and both are supported; the
caption must not let this panel be read as contradicting them.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import nm_style as S
import remote_data as D
from matplotlib.patches import Rectangle


def _rect(x, y, w, h, color, alpha):
    return Rectangle((x, y), w, h, facecolor=color, edgecolor="none",
                     alpha=alpha, zorder=0)

# ---- data ----------------------------------------------------------------
orc = D.oracle().set_index("scope")
# best-model PDS per gene, canonical multi-seed harness (D.best_model, committed CSV)
BEST_MODEL_PDS = {g: round(v, 2) for g, v in D.best_model().items()}

genes = S.GENE_ORDER  # TP53, KRAS, GATA1, JAK1
xs = {g: float(orc.loc[g, "PDS_oracle"]) for g in genes}
ys = {g: BEST_MODEL_PDS[g] for g in genes}
# 95% bootstrap interval on each ceiling, from the same committed CSV. This is what
# replaces the invented X_SEP: a ceiling whose interval covers 0.50 has not been
# shown to exceed chance, which is the paper's own measurement-limited criterion.
xlo = {g: float(orc.loc[g, "ci_lo"]) for g in genes}
xhi = {g: float(orc.loc[g, "ci_hi"]) for g in genes}

# ---- figure --------------------------------------------------------------
S.apply_rcparams()
fig, ax = S.panel(45.2, 43.9)

# x now reaches 0.42 so the lower bound of every ceiling interval is on-scale
# (KRAS 0.440 is the smallest); y spans the models plus the stretch of the y = x
# diagonal the shading needs.
X_LO, X_HI = 0.42, 0.86
Y_LO, Y_HI = 0.46, 0.72
CHANCE = 0.5

# The ONLY x threshold drawn is chance, which is defined throughout the paper.
BAND = 0.035  # half-width in PDS of the "model tracks its ceiling" band (in caption)

xg = np.linspace(X_LO, X_HI, 400)
lo_edge = np.clip(xg - BAND, Y_LO, Y_HI)
hi_edge = np.clip(xg + BAND, Y_LO, Y_HI)

# ceiling at or below chance: the benchmark cannot score a model here at all
ax.add_patch(_rect(X_LO, Y_LO, CHANCE - X_LO, Y_HI - Y_LO, S.LIGHT_GREY, 0.55))

# model tracks its ceiling (benchmark-solvable): a band hugging y = x
ax.fill_between(xg, lo_edge, hi_edge, facecolor=S.LIGHT_GREY, alpha=0.45,
                edgecolor="none", zorder=0)

# model below its ceiling (computation-limited), only where the ceiling clears chance
ax.fill_between(xg, Y_LO, lo_edge, where=(xg >= CHANCE), facecolor=S.LIGHT_GREY,
                alpha=0.18, edgecolor="none", zorder=0)

# y = x (model == ceiling). Solid hairline, direct-labelled: the (3,2) grey dash is
# reserved for chance across the whole paper and must not also mean "equality" here.
ax.plot([Y_LO, Y_HI], [Y_LO, Y_HI], color=S.GREY, lw=0.6, zorder=1)
# chance reference lines, house convention
ax.axvline(CHANCE, color=S.GREY, lw=0.6, ls=(0, (3, 2)), zorder=1)
ax.axhline(CHANCE, color=S.GREY, lw=0.6, ls=(0, (3, 2)), zorder=1)
ax.text(CHANCE - 0.011, Y_HI - 0.006, "chance", ha="right", va="top",
        fontsize=5.5, color=S.GREY, rotation=90, rotation_mode="anchor", zorder=2)

# Region labels name what is DRAWN; the caption maps them to the three regimes.
# The left strip (ceiling not above chance) is only 0.08 PDS wide and cannot hold a
# legible label at 5.5 pt, so it is named in the caption and identified on the page
# by the chance line it abuts.
ax.text(0.852, 0.468, "model below\nits ceiling", ha="right", va="bottom",
        fontsize=5.5, color=S.GREY, linespacing=1.1, zorder=2)
ax.text(0.648, 0.686, "model tracks\nits ceiling", ha="center", va="center",
        fontsize=5.5, color=S.GREY, linespacing=1.1, zorder=2)

# gene points, with the 95% bootstrap interval on the ceiling. Whether that
# interval crosses the chance line is the measurement-limited verdict.
for g in genes:
    ax.plot([xlo[g], xhi[g]], [ys[g], ys[g]], color=S.GENE_COLORS[g], lw=0.9,
            solid_capstyle="butt", zorder=4)
    ax.scatter(xs[g], ys[g], s=26, color=S.GENE_COLORS[g],
               edgecolor="white", linewidth=0.5, zorder=5, clip_on=False)

# Gene labels: centred above each point, except GATA1, which sits on the horizontal
# chance line and so is labelled below. Every point now carries a horizontal CI, so a
# side-anchored label would land on its own whisker or on a neighbour's.
# Label style is identical to panels f and g: 6 pt, gene colour, not bold.
label_off = {
    "TP53":  (0.0, 0.011, "bottom"),
    "KRAS":  (0.0, 0.011, "bottom"),
    "GATA1": (0.0, -0.011, "top"),
    "JAK1":  (0.0, 0.011, "bottom"),
}
for g in genes:
    dx, dy, va = label_off[g]
    ax.text(xs[g] + dx, ys[g] + dy, g, ha="center", va=va, fontsize=6,
            color=S.GENE_COLORS[g], zorder=6)

# axes
ax.set_xlim(X_LO, X_HI)
ax.set_ylim(Y_LO, Y_HI)
ax.set_xticks([0.5, 0.6, 0.7, 0.8])
ax.set_yticks([0.5, 0.6, 0.7])
ax.set_xlabel("Oracle ceiling (PDS, native depth) with 95% CI")
ax.set_ylabel("Best-model PDS")
S.despine(ax, keep=("left", "bottom"))

S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig5c_phasemap"))
print("done")
print("x (oracle):", {g: round(xs[g], 3) for g in genes})
print("y (best-model):", ys)
