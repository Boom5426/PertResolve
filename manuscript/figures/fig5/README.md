# Figure 5 — assembly guide

**One-line message:** the measurement window sets an achievable ceiling and governs
benchmark validity — a benchmark cannot rank models reliably when the ground-truth
discrimination ceiling is itself near chance. Source of truth:
`manuscript/latex/AllelePerturb_manuscript.tex`. House palette via shared
`manuscript/figures/nm_style.py`.

## Panel status (this pass)

| Panel | Content | Status | File / source |
|-------|---------|--------|---------------|
| a | oracle-ceiling definition (2nd measurement) | **done (prompt)** | `fig5a_prompt.md` (AI schematic) |
| d | synthetic-predictor construction (α interpolation) | **done (prompt)** | `fig5d_prompt.md` (AI schematic) |
| e | low- vs high-resolution leaderboard example | **done (prompt)** | `fig5e_prompt.md` (AI schematic) |
| f | **keystone**: benchmark resolution vs oracle ceiling | **done (data)** | `fig5f_resolution.py` — `benchmark_resolution/summary.csv` |
| b | per-gene oracle ceiling vs best model (dumbbell) | **HELD** | native-depth oracle 0.485/0.500/0.572/0.792 + CIs and best-model PDS are remote (Fig 2 bundle); local benchmark_resolution has a different ≥100-cell oracle (0.492/0.499/0.52/0.887, no CI) |
| c | measurement-limited vs computation-gap phase map | **HELD** | same native-oracle + best-model as 5b |
| g | P(select true best model) vs ceiling | **HELD** | committed files store only P(recover order), no winner-selection probability |

## 5f — verified numbers (self-contained, committed)

`benchmark_resolution/summary.csv`, (oracle_ceiling, resolution_P_recover_order),
reproduces SI Table 5 exactly:
TP53 0.492→0.014, KRAS 0.499→0.036, GATA1 0.520→0.512, sci-Plex 0.602→0.772,
VCC 0.731→1.00, Replogle 0.783→1.00, Adamson 0.836→1.00, Norman 0.848→1.00,
JAK1 0.887→0.99. Below a ceiling near 0.5 (TP53/KRAS) resolution collapses; above ~0.6
(GATA1 partial, then JAK1 + gene-/drug-level screens) it reaches ~1. The x-axis here is
the ≥100-cell / 50-per-half benchmark-resolution oracle (SI Table 5 depth), NOT the
native-depth per-gene oracle of the main-text window table (that is panel 5b, remote).

## Why 5b / 5c / 5g are held

- 5b/5c need the **native-depth per-gene oracle ceiling** (0.485/0.500/0.572/0.792) with
  95% CIs and the **best-model PDS** (0.51/0.52/0.54/0.52). The native oracle + CIs are
  not in any committed file, and best-model PDS is the same multi-seed/aggregation issue
  as Fig 2. Pull with the Fig 2 remote bundle (`unified/` + oracle native run).
- 5g needs a **P(select best model)** curve; committed files store only
  `resolution_P_recover_order` (full-order recovery), no winner-selection field.

The manuscript's regime-binned statements (ceiling≤0.52 → P=0.27; 0.52–0.56 → 0.95;
>0.65 → ~1.0) come from a finer controlled-predictor sweep (`controlled_recovery.csv`,
not committed); 5f here uses the committed per-dataset points, which tell the same story.

## Notes

- 5f imports the shared `nm_style`; vector PDF (0 raster) + 600 dpi PNG.
- AI panels (5a/5d/5e): paste prompt into a raster model, overlay exact text in
  Illustrator. Keep the no-results-leakage guard (no ceiling values / model names in the
  schematics).
- Do NOT put split-half ratio, external model audit, feature comparison, rankability
  predictor or the pilot workflow in Fig 5 (those are Fig 3/4/6).
