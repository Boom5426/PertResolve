"""Fig 5c — measurement-limited vs computation-gap phase map.

Message: two distinct failure modes separate in ceiling-vs-model space.
  x = oracle discrimination ceiling per gene (independent second measurement).
  y = best honestly-scored model PDS per gene.

Source of truth (both axes committed):
  x: results/_remote/unified/oracle_ceiling.csv via D.oracle() PDS_oracle
     (TP53 0.485, KRAS 0.500, GATA1 0.572, JAK1 0.792).
  y: manuscript Table tab:window "Best model PDS" column
     (AllelePerturb_manuscript.tex L226-229; best excluding PerturbNet, canonical
     multi-seed harness): TP53 0.51, KRAS 0.52, GATA1 0.54, JAK1 0.52.
The y column is committed only in the .tex table (not in any D-exposed CSV), so it
is pinned here verbatim with the trace above. No number is fabricated.

Regions (faint labels only, no leakage of ceiling/model values):
  lower-left  "measurement-limited"  both near chance (TP53, KRAS)
  lower-right "computation-limited"   ceiling high, model low (JAK1)
  upper-right "benchmark-solvable"    ceiling high, model tracks ceiling (aspiration)
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
# best-model PDS per gene, committed in manuscript Table tab:window (L226-229)
BEST_MODEL_PDS = {"TP53": 0.51, "KRAS": 0.52, "GATA1": 0.54, "JAK1": 0.52}

genes = S.GENE_ORDER  # TP53, KRAS, GATA1, JAK1
xs = {g: float(orc.loc[g, "PDS_oracle"]) for g in genes}
ys = {g: BEST_MODEL_PDS[g] for g in genes}

# ---- figure --------------------------------------------------------------
S.apply_rcparams()
fig, ax = S.panel(50.0, 48.0)

LO, HI = 0.45, 0.85
CHANCE = 0.5

# faint region shading: lower band (measurement-limited) uses ceiling <= ~0.55.
# We separate the plane by a vertical measurement threshold and the y=x diagonal.
X_SEP = 0.55  # ceiling below this => measurement-limited (near-chance measurement)

# lower-left measurement-limited block (ceiling near chance)
ax.add_patch(_rect(LO, LO, X_SEP - LO, HI - LO, S.LIGHT_GREY, 0.30))
# right half (high ceiling) split by the y=x diagonal into two failure/success modes:
#   below diagonal  -> computation-limited (ceiling high, model lags)
#   above/near diag -> benchmark-solvable  (model tracks ceiling)
comp_poly = [(X_SEP, LO), (HI, LO), (HI, HI), (X_SEP, X_SEP)]
solv_poly = [(X_SEP, X_SEP), (HI, HI), (X_SEP, HI)]
ax.add_patch(Polygon(comp_poly, closed=True, facecolor=S.GENE_COLORS["JAK1"],
                     edgecolor="none", alpha=0.12, zorder=0))
ax.add_patch(Polygon(solv_poly, closed=True, facecolor=S.GENE_COLORS["GATA1"],
                     edgecolor="none", alpha=0.10, zorder=0))

# y = x diagonal (model == ceiling)
ax.plot([LO, HI], [LO, HI], color=S.GREY, lw=0.6, ls=(0, (4, 2)), zorder=1)
# chance reference lines
ax.axhline(CHANCE, color=S.LIGHT_GREY, lw=0.5, zorder=0)
ax.axvline(CHANCE, color=S.LIGHT_GREY, lw=0.5, zorder=0)

# region labels (faint, no numbers)
ax.text(0.492, 0.472, "measurement-\nlimited", ha="center", va="center",
        fontsize=5.4, color=S.GREY, linespacing=1.0, zorder=2)
ax.text(0.775, 0.485, "computation-\nlimited", ha="center", va="center",
        fontsize=5.4, color=S.GREY, linespacing=1.0, zorder=2)
ax.text(0.635, 0.795, "benchmark-\nsolvable", ha="center", va="center",
        fontsize=5.4, color=S.GENE_COLORS["GATA1"], linespacing=1.0, zorder=2)

# gene points
for g in genes:
    ax.scatter(xs[g], ys[g], s=26, color=S.GENE_COLORS[g],
               edgecolor="white", linewidth=0.5, zorder=5, clip_on=False)

# gene labels, offset to avoid overlap and stay on-axis
# TP53 (0.485,0.51) & KRAS (0.500,0.52) sit close in lower-left -> spread them
label_off = {
    "TP53":  (-0.006, 0.018, "right", "bottom"),
    "KRAS":  (0.010, -0.006, "left", "top"),
    "GATA1": (0.012, 0.004, "left", "center"),
    "JAK1":  (-0.010, 0.006, "right", "bottom"),
}
for g in genes:
    dx, dy, ha, va = label_off[g]
    ax.text(xs[g] + dx, ys[g] + dy, g, ha=ha, va=va, fontsize=6,
            color=S.GENE_COLORS[g], fontweight="bold", zorder=6)

# axes
ax.set_xlim(LO, HI)
ax.set_ylim(LO, HI)
ax.set_xticks([0.5, 0.6, 0.7, 0.8])
ax.set_yticks([0.5, 0.6, 0.7, 0.8])
ax.set_xlabel("Oracle ceiling (PDS)")
ax.set_ylabel("Best-model PDS")
ax.set_aspect("equal", adjustable="box")
S.despine(ax, keep=("left", "bottom"))

S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig5c_phasemap"))
print("done")
print("x (oracle):", {g: round(xs[g], 3) for g in genes})
print("y (best-model):", ys)
