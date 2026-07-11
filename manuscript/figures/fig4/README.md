# Figure 4 — assembly guide

**One-line message:** allele-ranking failure persists across evaluation choices and
model interfaces (splits, metrics, representations, model class) — an
alternative-explanation-elimination figure. Source of truth:
`manuscript/latex/AllelePerturb_manuscript.tex`. House palette via shared
`manuscript/figures/nm_style.py`.

## Panel status (this pass)

| Panel | Content | Status | File / source |
|-------|---------|--------|---------------|
| c | empirical permutation-null calibration | **done** | `fig4c_null.py` — `permutation_null_pds.csv` (self-contained; verified null mean 0.500, 3.2% exceed) |
| h | interface-compatibility audit matrix | **done** | `fig4h_interface.py` — categorical, from `docs/REMOTE_INVENTORY.md` + Methods |
| a | PDS across 5 splits | **HELD** | needs the canonical multi-seed PDS (remote `unified/`, same source as Fig 2) |
| b | Pearson-δ across 5 splits | **HELD** | pairs with 4a; Pearson from grid, PDS from `unified/` |
| d | distance invariance (cosine/L1/L2) | **HELD** | multi-seed PDS per distance (remote `unified/`) |
| e | feature spaces (θ / ESM / ESM+θ) | **HELD** | multi-seed PDS + Pearson (remote `unified/`) |
| f | gene × split dissociation | **HELD** | multi-seed PDS per gene×split (remote `unified/`) |
| g | **hero**: published external models under one harness | **HELD** | no local external-model results (remote `unified/preds5`, `definitive_summary.csv`) |

**Why 4a/b/d/e/f/g are held:** their PDS values must come from the SAME canonical
multi-seed scorer chosen for Fig 2 (remote `unified/`), so all figures quote one
scorer. The local single-draw grid (`results_v4_10metrics.csv`, what the legacy
`draw_fig4.py` uses) gives PDS ~0.44–0.50, not the manuscript's 0.49–0.52. Server was
down 2026-07-11; pull with the Fig 2 bundle when it returns.

## 4c — verified numbers (self-contained, no remote dependency)

`permutation_null_pds.csv`, metric PDS_cos, n = 340 method×split×gene combinations:
permutation null mean **0.500** (per-combo means 0.495–0.504); **3.2%** exceed their
own 95th-percentile null (expected 5%); **33%** of combinations land exactly at 0.50
(tie-dominated small test sets). This closes the "maybe 0.50 is the wrong chance
baseline" objection.

## 4h — interface audit content (categorical facts)

Rows grouped as: **variant-conditionable** (scGen, scVIDR, Biolord, CellFlow ran to
completion and are scored; PerturbNet ran but its score is subspace-caveated) ·
**allele-blind by construction** (scGPT, GEARS, STATE condition on gene identity;
STATE was run per-gene but is still allele-blind) · **did not converge** (variant-CPA:
continuous θ input attempted, non-finite loss). Columns: continuous variant input /
unseen-allele conditioning / allele-specific output / ran to completion / allele-level
score defined. Gene-keyed resolution coverage = **4/470** variants (0.85%) — updated
from 4/472 with the WT-row variant-count correction. Framing note: gene-keyed models
are an *interface incompatibility* (allele-level score undefined), NOT a performance
failure — keep that wording.

## Notes

- 4c and 4h import the shared `nm_style`; each writes a vector PDF (0 embedded raster)
  + 600 dpi PNG. 4h glyphs (check/cross/tilde) are drawn as shapes (Liberation Sans
  lacks ✓/✗), so they are font-independent.
- The hero panel 4g and the robustness panels are the bulk of Fig 4 and depend on the
  Fig 2 remote sync; only calibration (4c) and the interface audit (4h) are drawable now.
