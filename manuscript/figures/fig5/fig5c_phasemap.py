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
  left band   "measurement-limited"  ceiling near chance (TP53, KRAS)
  lower-right "computation-limited"   ceiling high, model low (JAK1)
  diagonal    "benchmark-solvable"    ceiling high, model TRACKS the ceiling, i.e. a
              narrow band hugging y = x in the high-ceiling zone. The half-plane
              ABOVE y = x (model scores higher than its own oracle ceiling) is not
              an attainable regime: it is left unshaded and unlabelled, since any
              excursion into it is best-of-many-heads selection noise (see panel b).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D
from matplotlib.patches import Polygon, Rectangle


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

# ---- figure --------------------------------------------------------------
S.apply_rcparams()
fig, ax = S.panel(45.2, 43.9)

# y-axis is tightened to the plotted models plus the stretch of the y = x diagonal
# that the regime shading needs; x still spans every gene ceiling (JAK1 0.792).
X_LO, X_HI = 0.46, 0.82
Y_LO, Y_HI = 0.46, 0.72
CHANCE = 0.5

# faint region shading. The plane is split by a vertical measurement threshold and
# by the y = x diagonal (model == its own ceiling).
X_SEP = 0.55  # ceiling below this => measurement-limited (near-chance measurement)
BAND = 0.035  # half-width (in PDS) of the "model tracks its ceiling" band

# left measurement-limited block (ceiling near chance)
ax.add_patch(_rect(X_LO, Y_LO, X_SEP - X_LO, Y_HI - Y_LO, S.LIGHT_GREY, 0.55))

# benchmark-solvable: a narrow band hugging y = x once the ceiling is high, i.e.
# the model tracks its ceiling. Clipped at the top of the axes.
x_top = min(Y_HI + BAND, X_HI)      # y = x - BAND reaches the top edge here
solv_poly = [(X_SEP, X_SEP - BAND), (x_top, Y_HI), (Y_HI - BAND, Y_HI),
             (X_SEP, X_SEP + BAND)]
ax.add_patch(Polygon(solv_poly, closed=True, facecolor=S.LIGHT_GREY,
                     edgecolor="none", alpha=0.45, zorder=0))

# computation-limited: high ceiling, model below the tracking band.
comp_poly = [(X_SEP, Y_LO), (X_HI, Y_LO), (X_HI, Y_HI), (x_top, Y_HI),
             (X_SEP, X_SEP - BAND)]
ax.add_patch(Polygon(comp_poly, closed=True, facecolor=S.LIGHT_GREY,
                     edgecolor="none", alpha=0.18, zorder=0))

# y = x diagonal (model == ceiling)
ax.plot([X_LO, Y_HI], [X_LO, Y_HI], color=S.GREY, lw=0.6, ls=(0, (4, 2)), zorder=1)
# chance reference lines
ax.axhline(CHANCE, color=S.LIGHT_GREY, lw=0.5, zorder=0)
ax.axvline(CHANCE, color=S.LIGHT_GREY, lw=0.5, zorder=0)

# region labels (faint, no numbers). "benchmark-solvable" runs ALONG the diagonal
# so it reads as "model tracks its ceiling", not as the half-plane above it.
ax.text(0.505, 0.680, "measurement-\nlimited", ha="center", va="center",
        fontsize=5.5, color=S.GREY, linespacing=1.0, zorder=2)
ax.text(0.660, 0.558, "computation-\nlimited", ha="center", va="center",
        fontsize=5.5, color=S.GREY, linespacing=1.0, zorder=2)
# centred ON the diagonal (x = y = 0.665), i.e. inside the tracking band, so the
# label names the band and not the half-plane above it. Kept horizontal: rotated
# text is fragmented by the on-page type audit (pdftotext) into sub-word boxes.
ax.text(0.665, 0.665, "benchmark-\nsolvable", ha="center", va="center",
        fontsize=5.5, color=S.GREY, linespacing=1.0, zorder=2)

# gene points
for g in genes:
    ax.scatter(xs[g], ys[g], s=26, color=S.GENE_COLORS[g],
               edgecolor="white", linewidth=0.5, zorder=5, clip_on=False)

# gene labels, offset to avoid overlap and stay on-axis
# TP53 (0.485,0.51) & KRAS (0.500,0.54) sit close in lower-left -> spread them
# label style is identical to panels f and g: 6 pt, gene colour, not bold.
# (dy rescaled by the tightened y-range so the on-page offsets are unchanged)
label_off = {
    "TP53":  (-0.007, 0.001, "right", "center"),
    "KRAS":  (0.009, 0.010, "left", "bottom"),
    "GATA1": (0.010, -0.003, "left", "top"),
    "JAK1":  (-0.008, 0.010, "right", "bottom"),
}
for g in genes:
    dx, dy, ha, va = label_off[g]
    ax.text(xs[g] + dx, ys[g] + dy, g, ha=ha, va=va, fontsize=6,
            color=S.GENE_COLORS[g], zorder=6)

# axes
ax.set_xlim(X_LO, X_HI)
ax.set_ylim(Y_LO, Y_HI)
ax.set_xticks([0.5, 0.6, 0.7, 0.8])
ax.set_yticks([0.5, 0.6, 0.7])
ax.set_xlabel("Oracle ceiling (PDS, native depth)")
ax.set_ylabel("Best-model PDS")
S.despine(ax, keep=("left", "bottom"))

S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig5c_phasemap"))
print("done")
print("x (oracle):", {g: round(xs[g], 3) for g in genes})
print("y (best-model):", ys)
