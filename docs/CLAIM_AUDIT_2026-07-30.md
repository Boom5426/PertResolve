# Claim audit, 2026-07-30

Step 4 of the Nature Methods Article conversion. Every empirical claim in the abstract,
introduction, title, Results subsections, figure legends, Discussion generalisations,
Methods reproducibility statements and Supplementary tables was traced to a committed
result file and recomputed where the value is derived.

**Status.** The audit itself edited nothing. A follow-up step-3 pass then applied the ten text
corrections below and left the AlphaMissense items untouched by decision:

| Item | State |
|---|---|
| B1 residue-cluster sentence | FIXED |
| W1 27% / 38% bin scope | FIXED |
| W2 GATA1 ceiling label (Results **and** the Fig. 5g legend, which carried the same mislabel) | FIXED |
| W3 the two 15-of-17 exceptions | FIXED |
| W4 Fig. 3g PerturbNet legend | FIXED |
| W5 Fig. 2g PDS range and fraction | FIXED |
| W8 titration space and JAK1 n (Results **and** the Fig. 4g legend) | FIXED |
| item 7 atlas duplication | FIXED, compressed in Results with the CIs moved to the pilot subsection |
| W9 JAK1 panel size in limitations | FIXED |
| B4 AlphaMissense artifact | **RESOLVED**: `am_summary.csv` and `am_effect_spearman.csv` committed, script de-hardcoded and now writes its tables |
| W6 Pearson-delta upper bound 0.66 | **PASS**: confirmed as AlphaMissense, and the range still rounds to 0.63 to 0.66 after recomputation |
| B4-bis MLP head silently omitted | **RESOLVED AFTER RECOMPUTATION**: found while fixing B4, see below; representation grid re-run on six heads, four text numbers updated |
| B2, B3, B5, pilot per-dataset JSON | open, code and repository work |

Both documents rebuild clean after the pass (0 errors, 0 undefined citations or references,
0 overfull boxes, 0 float warnings; main 26 pages, SI 7 pages).

**B4 follow-up (2026-07-30, later the same day).** Producing the AlphaMissense artifact
uncovered a defect that the bare `except Exception` had been hiding. See B4-bis below. The
AlphaMissense values are now committed and confirm the reported ranges, but the representation
grid had to be re-run, and four numbers need updating as a result.

## Evidence base

| Layer | Location | Status |
|---|---|---|
| Benchmark composition | `data/allele_perturb_bench.csv` | tracked |
| Canonical multi-seed tables | `results/canonical/*.csv`, `*.json` | tracked, mapped by `results/canonical/README.md` |
| Single-draw grid | `results/results_v4_10metrics.csv`, `results/results_v4_exttheta.csv` | tracked |
| Reviewer controls | `results/reviewer_controls/` | tracked (tables + scripts) |
| Pilot validation | `results/pilot_validation/` | tracked (tables + scripts + README) |
| Benchmark resolution | `results/benchmark_resolution/` | tracked (tables + scripts) |
| Local mirror of remote outputs | `results/_remote/` | gitignored by design; canonical subset promoted to `results/canonical/` |

Aggregation conventions are not interchangeable. `manuscript/figures/remote_data.py`
documents a dual source and the audit used it:

* per-method forests and the external-model table: `results/canonical/definitive_summary.csv`
  and `unified_multiseed.csv` (multi-seed, 15 subsamples);
* per-gene / per-split / per-variant / per-feature breakdowns: `results/results_v4_exttheta.csv`.

An earlier pass of this audit used `unified_results5.csv` for the breakdown panels and
produced four spurious mismatches. The panel scripts in `manuscript/figures/fig2/` are the
operative definition and were used instead.

## Verdict summary

| Verdict | n | Severity |
|---|---|---|
| PASS | 45 | |
| ARTIFACT_MISSING | 5 | submission blocker |
| VALUE_MISMATCH | 4 | must fix before submission |
| WORDING_TOO_STRONG | 5 | must fix before submission |

No claim was found where the underlying science is wrong. Every defect is a
traceability gap or a scope/label imprecision.

---

# Blockers

## B1. Methods residue-cluster claim is contradicted by its own committed table

```text
Location:        Methods, Statistical analysis
Exact claim:     "we also resampled residue-level clusters rather than individual variants.
                 The per-method intervals were essentially unchanged, with every non-null
                 predictor still overlapping chance."
Claim type:      methodological / numerical
Supporting file: results/reviewer_controls/residue_independence_summary.csv
Exact row:       PerturbNet, split1: pds_mean 0.565, pds_clust_lo 0.502, pds_clust_hi 0.625,
                 pds_overlaps_0p5 = False
Verified value:  23 of 24 methods overlap 0.50. PerturbNet does not (cluster CI 0.502-0.625).
                 The table also covers split1 only (24 rows, one split), so "the per-method
                 intervals" is a random-split statement, not a global one.
Verdict:         VALUE_MISMATCH
Required action: carve out PerturbNet and name the split. This does not weaken the paper:
                 PerturbNet is already the single non-crossing method elsewhere
                 (definitive_summary.csv crosses=False) and is already explained as a
                 subspace artifact. Suggested wording:
                 "Resampling residue-level clusters rather than individual variants left the
                 random-split intervals essentially unchanged: 23 of 24 predictors still
                 overlapped chance, the exception being PerturbNet, whose exceedance is the
                 50-dimensional subspace artifact characterised above."
```

## B2. Code availability overstates what the public repo regenerates

```text
Location:        Methods, Code availability
Exact claim:     "All analysis code, the packaged benchmark (alleleperturb/) and the figure
                 and evaluation scripts that regenerate the main-figure and Supplementary
                 values are available at https://github.com/Boom5426/AllelePerturb."
Claim type:      methodological
Supporting file: git ls-files scripts/ alleleperturb/ ; results/canonical/README.md
Verified value:  9 of 11 canonical tables have no committed generating script:
                 definitive_summary.csv, unified_multiseed.csv, controlled_recovery.csv,
                 pairwise_resolvability.csv, classifier_two_sample.csv,
                 metric_family_floor.csv, floor_law.csv, leakage_control_summary.csv,
                 best_model_pds.csv.
                 The figure scripts regenerate the panels from these CSVs; the harness that
                 computed the CSVs is not in the repo.
Verdict:         ARTIFACT_MISSING
Required action: either commit the canonical harness, or narrow the sentence to what is true
                 ("the figure scripts regenerate every main-figure panel from the committed
                 canonical tables in results/canonical/") and release the harness alongside
                 the packaging work. Folded into the software fix list below.
```

## B3. PerturbNet permutation P = 0.13 has no committed generating code

```text
Location:        Results (metrics/representations/interfaces subsection);
                 Fig. 3g legend; Supplementary Table 4 footnote
Exact claim:     "against a permutation null estimated in that subspace its exceedance was
                 not significant (permutation P = 0.13)"
Claim type:      numerical
Supporting file: results/canonical/perturbnet_subspace_test.json
Exact field:     perturbnet_permutation_p = 0.129; perturbnet_permutation_null_mean = 0.5435;
                 perturbnet_subspace_pds = 0.561; pure_random_in_subspace_pds = 0.516
Verified value:  0.129 -> 0.13. Value present in a tracked artifact.
                 BUT scripts/analysis/subspace_test.py contains no permutation test at all
                 (only the four projection scores), writes nothing (print only, no to_csv /
                 json.dump), and hard-codes BASE = "/data/boom/NUS/VCCompass".
Verdict:         ARTIFACT_MISSING (code)
Required action: extend subspace_test.py to compute and write the permutation test with an
                 explicit seed and stated n_permutations; de-hardcode the data path.
```

## B4. AlphaMissense is claimed as a scored representation with no committed result

```text
Location:        Results (metrics/representations/interfaces subsection); SI representation list
Exact claim:     "...and AlphaMissense, a purpose-built structure- and evolution-informed
                 variant-effect predictor. Every representation left allele discrimination at
                 chance (PDS 0.46 to 0.48)"
                 and "AlphaMissense pathogenicity did not track the measured transcriptional
                 effect size within any gene (Spearman near zero for TP53, KRAS and GATA1)"
Claim type:      numerical / scope
Supporting file: results/reviewer_controls/rep_grid_summary.csv (7 rows: theta, esm1v,
                 esm2_global, esm2_globaldelta, esm2_sitedelta, esm2_window16,
                 esm2_sitedelta+theta) -- no AlphaMissense row.
                 results/reviewer_controls/am_analysis.py DOES compute both quantities
                 (section 1 Spearman, section 2 "AM as a feature through the 6 heads") but
                 only prints them; no to_csv. It hard-codes BASE = "/data/boom/NUS/VCCompass"
                 and reads AM_4genes.tsv, which is not in the repo.
Verified value:  the 0.46-0.48 PDS range is exact for the 7 committed representations
                 (0.464-0.475). AlphaMissense is inside the claim's scope but outside its
                 evidence.
Verdict:         ARTIFACT_MISSING
Required action: write am_analysis.py outputs to results/reviewer_controls/am_summary.csv
                 (per-gene Spearman + AM-as-feature PDS and Pearson-delta), or remove
                 AlphaMissense from the "every representation ... PDS 0.46 to 0.48" range and
                 report it only as the pathogenicity-versus-effect-size correlation.
Note:            this is also the root cause of W6. The Pearson-delta upper bound 0.66 is
                 AlphaMissense's 0.656; the 2026-07-12 pass raised 0.65 to 0.66 for exactly
                 that reason. Committing the AlphaMissense row fixes both entries at once.
RESOLVED:        2026-07-30. am_analysis.py now takes --base / VCCOMPASS_BASE and writes
                 results/reviewer_controls/am_summary.csv (alphamissense 0.458 / 0.659,
                 alphamissense+theta 0.465 / 0.657) and am_effect_spearman.csv (per-gene
                 Spearman -0.0752, -0.0538, 0.0874, 0.4091 plus a provenance-only POOLED row).
                 Both claims are now traceable. Fixing this uncovered B4-bis below.
```

## B4-bis. The representation control ran on five heads, not six

Found while producing the B4 artifact. Recorded as a separate entry rather than folded into B4,
because it is an additional defect with its own scope, not a restatement of the artifact gap.

```text
Location:        Results, representations sentence ("Under the identical splits, heads and
                 scoring"); Methods, Regression heads ("six standard regression heads ... and a
                 multilayer perceptron"); Supplementary Note, Representation robustness
Claim type:      methodological
Root cause:      results/reviewer_controls/rep_grid.py and am_analysis.py both constructed the
                 MLP as MLPRegressor((128, 64), max_iter=500, random_state=0). In
                 scikit-learn >= 1.7 the first POSITIONAL parameter of MLPRegressor is `loss`,
                 not `hidden_layer_sizes`, so every fit raised
                   InvalidParameterError: The 'loss' parameter of MLPRegressor must be a str
                   among {'poisson','squared_error'}. Got (128, 64) instead.
                 and a bare `except Exception: continue` discarded it. The failure rate was
                 16 of 16 evaluated gene-by-split cells, i.e. the neural-network family
                 contributed nothing to any of the seven representations or to AlphaMissense.
                 Confirmed on the compute host: scikit-learn 1.7.2, keyword form fits, the
                 positional form raises.
Verified value:  re-run 2026-07-30, six heads, zero failures. Every value moved by at most
                 0.003 and no verdict changed.

                   representation          PDS_cos 5h -> 6h    Pearson-delta 5h -> 6h
                   theta                   0.475 -> 0.473      0.632 -> 0.632
                   esm1v                   0.464 -> 0.462      0.640 -> 0.643
                   esm2_global             0.465 -> 0.463      0.645 -> 0.648
                   esm2_globaldelta        0.465 -> 0.463      0.646 -> 0.648
                   esm2_sitedelta          0.468 -> 0.467      0.631 -> 0.632
                   esm2_window16           0.466 -> 0.464      0.646 -> 0.648
                   esm2_sitedelta+theta    0.465 -> 0.464      0.627 -> 0.628
                   alphamissense           0.459 -> 0.458      0.655 -> 0.659
                   alphamissense+theta     0.466 -> 0.465      0.656 -> 0.657

                 Nine-representation range: PDS-cosine 0.458 to 0.473, Pearson-delta 0.628 to
                 0.659. Every representation remains at chance.
Scope:           the MAIN evaluation grid is NOT affected. scripts/run_all_splits.py:181-183
                 and scripts/figures/all_splits_v4.py:181-183 use the keyword form, and the
                 committed MLP-theta, MLP-esm and MLP-esm+theta rows carry real values
                 (0.495, 0.495, 0.493). Fig. 2, Fig. 3, the 18 in-house heads, every headline
                 PDS and Pearson-delta, the oracle ceilings and the benchmark-resolution
                 analysis are unchanged.
Verdict:         RESOLVED AFTER RECOMPUTATION
Code fixes:      both scripts pass hidden_layer_sizes by keyword with a comment naming the
                 cause; both take --base / VCCOMPASS_BASE instead of a hard-coded
                 /data/boom/NUS/VCCompass; both write their tables to --out; both count
                 head-fit failures per gene, split and head instead of swallowing them.
                 rep_grid_summary.csv gained n_heads_scored, heads and n_head_failures so the
                 omission cannot recur silently.
Provenance:      the five-head table is archived at
                 results/deprecated/rep_grid_summary_5heads.csv with the full before-and-after
                 comparison; results/reviewer_controls/README.md carries the caveat.
Text updated:    Results "PDS 0.46 to 0.48" -> "0.46 to 0.47" (Pearson-delta 0.63 to 0.66
                 unchanged, since 0.628 to 0.659 still rounds to it); Supplementary Note
                 "0.459--0.475" -> "0.458--0.473", "0.627--0.656" -> "0.628--0.659",
                 gene-mean "at PDS 0.45" -> "at PDS 0.44" (results_v4_10metrics.csv 0.4416,
                 results_v4_exttheta.csv 0.4432, both round to 0.44), sources extended to
                 am_summary.csv and am_effect_spearman.csv, and a coverage clause added.
Coverage:        the scoring procedure is identical across representations but the variant
                 coverage is not, because each representation exists only where its input
                 does. AlphaMissense covers 374 of the 470 benchmark variants (TP53 84/98,
                 KRAS 78/92, GATA1 201/254, JAK1 11/26), i.e. 389 scored variant-by-split
                 instances against 518 for theta, 511 for ESM-1v and 497 for the ESM2
                 constructions. Now stated in the Supplementary Note.
KRAS universe:   unified/real_deltas.npz holds 97 KRAS keys, five more than the curated
                 benchmark's 92 (AG11TD, AG59GV, C185Y, K179R, M170L). They are absent from
                 allele_perturb_bench.csv, which is what harness.split_vars reads, so they
                 never enter a train, test or candidate set, and none carries an AlphaMissense
                 score. Every reported denominator is therefore the benchmark universe.
Pooling red line: am_effect_spearman.csv keeps a POOLED row (-0.3798) that disagrees in sign
                 with three of the four per-gene values (-0.0752, -0.0538, 0.0874, 0.4091).
                 Cross-gene pooling is not interpretable, because both the pathogenicity
                 distribution and the measured effect scale differ between genes, so the
                 pooled value is driven by between-gene offsets. DO NOT REPORT it. JAK1's
                 0.409 rests on 11 of 26 variants and is likewise not a positive result. The
                 manuscript reports the within-gene TP53, KRAS and GATA1 values only.
```

## B5. The released package does not implement the framework the abstract advertises

```text
Location:        Title; abstract; Methods Code availability
Exact claim:     "a power-aware framework that separates perturbation detection,
                 sibling-allele identification and held-out response prediction, estimates
                 replicate-based discrimination ceilings, and predicts rankability from pilot
                 measurements"
Claim type:      methodological
Supporting file: alleleperturb/ (495 lines: bench.py, features.py, metrics.py,
                 evaluation/{pds,direction_metrics,de_metrics}.py)
Verified value:  the package implements benchmark loading, splits, theta, and three metric
                 families. Split-half detection diagnostics (D_self/D_null), sibling-allele
                 identification, the replicate oracle ceiling and pilot rankability
                 prediction are absent. No CLI, no def main, no entry_points, no pyproject.toml.
Verdict:         ARTIFACT_MISSING
Required action: lift the four diagnostics into the package (diagnostics.py, ceiling.py,
                 rankability.py) with a CLI and pyproject.toml. This is the standing
                 Nature Methods Article blocker, not a prose problem.
```

---

# Wording and scope corrections

## W1. The 27% / 38% bin is mislabelled

```text
Location:        Results, replicate-ceilings subsection
Exact claim:     "When the oracle window was at chance (ceiling <= 0.52; TP53 and most depths
                 of KRAS), the benchmark recovered the true predictor order only 27% of the
                 time and identified the best predictor only 38% of the time"
Claim type:      numerical / scope
Supporting file: results/canonical/controlled_recovery.csv
Exact rows:      the 7 rows with ceiling_pds <= 0.52 are TP53@25, TP53@50, KRAS@25, KRAS@50,
                 KRAS@100, KRAS@150 and GATA1@25 (ceiling 0.517, P_correct 0.872)
Verified value:  mean over those 7 rows = P_correct 0.2713, P_winner 0.3839 -> 27% / 38%,
                 exact. TP53 + KRAS only = 17% / 30%.
Verdict:         WORDING_TOO_STRONG (scope label)
Required action: "(ceiling <= 0.52, comprising TP53, all KRAS depths and the shallowest GATA1
                 configuration)". Note the correction is conservative: restricting the bin to
                 TP53 and KRAS lowers the numbers to 17% and 30%.
```

## W2. "ceiling ~0.57; GATA1" labels the wrong quantity

```text
Location:        Results, replicate-ceilings subsection
Exact claim:     "Once the window opened even modestly (ceiling ~0.57; GATA1), correct
                 ordering rose to 95% and correct winner selection to 97%"
Claim type:      numerical
Supporting file: results/canonical/controlled_recovery.csv; results/canonical/oracle_ceiling.csv
Verified value:  0.9455 and 0.9692 reproduce exactly as the mean over GATA1 depths
                 50/100/150/250, whose ceilings are 0.522 to 0.549. The value 0.572 is
                 GATA1's native-depth oracle ceiling from oracle_ceiling.csv, a different
                 quantity from the configurations being averaged.
Verdict:         VALUE_MISMATCH (label)
Required action: "(ceiling 0.52 to 0.55; GATA1)".
```

## W3. "both metrics were close" is false for one of the two exceptions

```text
Location:        Results, metrics/representations/interfaces subsection
Exact claim:     "Pearson-delta exceeded PDS in 15 of 17 combinations, with the two
                 exceptions occurring in small held-out mechanistic-extrapolation sets where
                 both metrics were close."
Claim type:      numerical / qualitative
Supporting file: results/results_v4_exttheta.csv, 18 in-house heads
Verified value:  exactly 17 gene-by-split combinations, 15 with Pearson > PDS. The two
                 exceptions are KRAS/split3 (Pearson 0.709 vs PDS 0.715, gap -0.006, close)
                 and JAK1/split3 (Pearson 0.022 vs PDS 0.393, gap -0.371, not close;
                 JAK1's direction recovery itself collapses).
Verdict:         WORDING_TOO_STRONG
Required action: "...the two exceptions occurring in small held-out mechanistic-extrapolation
                 sets: KRAS, where the two metrics coincide (0.71 versus 0.72), and JAK1,
                 where direction recovery itself collapses (Pearson-delta 0.02)."
```

## W4. The Fig. 3g legend asserts a containment the numbers do not show

```text
Location:        Fig. 3g legend (formerly Fig. 4g)
Exact claim:     "PerturbNet's point estimate (0.54) exceeds chance in the full harness but
                 falls inside its own 50-dimensional subspace null (0.52)"
Claim type:      numerical
Supporting file: results/canonical/perturbnet_subspace_test.json
Verified value:  0.54 > 0.52, so it does not fall inside 0.52. The precise statement is the
                 one already in the Supplementary Table 4 footnote: random predictions
                 confined to the subspace already reach 0.52, and PerturbNet's in-subspace
                 exceedance over a subspace permutation null is not significant
                 (in-subspace PDS 0.561, null mean 0.543, P = 0.13).
Verdict:         WORDING_TOO_STRONG
Required action: align the legend to the SI footnote.
```

## W5. Fig. 2g per-variant PDS upper bound and "roughly half"

```text
Location:        Results, direction-versus-discrimination subsection; Fig. 2g legend
Exact claim:     "PDS varied from 0.00 to 0.98, with roughly half of examples below chance"
Claim type:      numerical
Supporting file: results/results_v4_exttheta.csv via the exact selection in
                 manuscript/figures/fig2/fig2g_pervariant.py (Ridge-esm, TP53, dedup by
                 variant after sorting on PDS_cos, 10 evenly spaced)
Verified value:  Pearson-delta 0.735 to 0.821 -> "0.74 to 0.82" PASS.
                 PDS 0.000 to 0.948 -> 0.95, not 0.98. 6 of 10 below chance.
Verdict:         VALUE_MISMATCH (minor)
Required action: "0.00 to 0.95" and "more than half".
```

## W6. Representation Pearson-delta upper bound depends on the uncommitted AlphaMissense row

```text
Location:        Results, metrics/representations/interfaces subsection
Exact claim:     "similar direction recovery (Pearson-delta 0.63 to 0.66)"
Claim type:      numerical
Supporting file: results/reviewer_controls/rep_grid_summary.csv
Verified value:  the seven committed representations span pearson_delta 0.627 to 0.646, i.e.
                 "0.63 to 0.65". The 0.66 upper bound is AlphaMissense's 0.656, which was
                 deliberately introduced in the 2026-07-12 pass (0.65 -> 0.66) and is correct
                 as a value, but AlphaMissense has no committed row: it is computed only by
                 the print-only section 2 of am_analysis.py (see B4).
                 The companion PDS claim "0.46 to 0.48" is exact for the seven committed
                 representations (0.464 to 0.475).
Verdict:         PASS (resolved with B4; never an arithmetic error)
RESOLVED:        2026-07-30. With am_summary.csv committed and the grid re-run on six heads the
                 nine-representation range is Pearson-delta 0.628 to 0.659, which still rounds
                 to "0.63 to 0.66", so that number stands unchanged. The companion PDS range
                 became 0.458 to 0.473 and its printed upper bound was narrowed from 0.48 to
                 0.47. See B4-bis.
```

## W7. One sentence mixes the two documented conventions

```text
Location:        Results, metrics/representations/interfaces subsection
Exact claim:     "a featureless predictor that outputs the training gene mean already reached
                 Pearson-delta 0.66"
Claim type:      numerical
Supporting file: results/canonical/unified_results5.csv (Gene-mean pearson_delta = 0.664),
                 matching Supplementary Table 4's Gene-mean 0.66
Verified value:  0.664 -> 0.66 PASS. But the surrounding numbers in the same paragraph come
                 from the single-draw representation grid, and the same quantity in
                 results_v4_exttheta.csv (the grid backing the neighbouring panels) is 0.642.
Verdict:         PASS with a source-precision caveat
Required action: name the harness once, or use 0.64 for internal consistency with the grid
                 the rest of the paragraph cites.
```

## W8. Depth-titration claim is space-dependent and rests on n = 1 at the top depth

```text
Location:        Results, split-half/sibling-diagnostics subsection; Fig. 4g legend
Exact claim:     "Subsampling over 50 to 300 cells per half confirmed the same ordering:
                 JAK1 remained fully detectable across the range, whereas TP53, KRAS and
                 GATA1 increased only modestly and plateaued below 75%."
Claim type:      numerical / scope
Supporting file: results/split_half_power_curve.csv
Exact rows:      the file stores two representation spaces. In pca50 (the canonical space
                 that manuscript/figures/fig3/fig3e_titration.py reads) the maxima are
                 TP53 0.690, KRAS 0.735, GATA1 0.573, JAK1 1.000, so "below 75%" holds.
                 In the `full` space TP53 reaches 0.897 and KRAS 0.898, so the claim fails.
                 The manuscript never names the space.
                 JAK1 n_variants at n_sub 50/100/150/300 = 14 / 9 / 6 / 1. TP53 has no
                 300-cell row at all.
Verdict:         WORDING_TOO_STRONG (scope)
Required action: state PCA-50 in the sentence or the legend, and give JAK1's per-depth n so
                 "fully detectable at 300 cells per half" is not read as a population claim
                 when it rests on one variant.
```

## W9. Limitations omit JAK1's variant count

```text
Location:        Discussion, limitations paragraph
Claim type:      scope
Verified value:  the model-limited verdict rests entirely on JAK1: 26 variants total, 20 in
                 the rankability analysis, 14 at matched depth, 6 to 9 at the higher titration
                 depths; the 0.792 oracle ceiling has n_var_scored = 59 variant-split
                 instances.
Verdict:         WORDING (missing scope)
Required action: add JAK1's n to the limitations paragraph. The paragraph already carries the
                 gene/assay confound, the pilot's detection-level restriction, the
                 oracle-is-not-an-upper-bound caveat and the pseudobulk caveat.
```

---

# The user's seven priority items

## 1. Extended Data Figure 1: RESOLVED, not a gap

`AllelePerturb_SI.tex:40` contains `\includegraphics[width=\textwidth]{figures/ED_fig1.pdf}`,
followed at line 41 by a caption ("The rankability verdict is robust to the choice of
criterion"). `manuscript/latex/figures/ED_fig1.pdf` is tracked in git. It renders as the
SI's Figure 1, which the main text cites twice as "Supplementary Fig. 1"; both citations
match the caption's content ("alternative rankability criteria, including bootstrap-CI-based,
increased-stringency and ratio-threshold definitions"). Numbers trace to
`results/rankability_sensitivity.csv`.

No orphaned reference, no text-claims-it-but-PDF-lacks-it state. The only residual is
cosmetic: the file is named `ED_fig1.pdf` but is numbered as a Supplementary Figure.

## 2. PerturbNet 0.52 subspace null: value traceable, code not

Not a hard-coded constant in prose. See B3. The projection controls are in a committed
script; the permutation test is not, and the script cannot run from the public repo. The SI
Table 4 footnote states the result precisely; the Fig. 3g legend does not (W4).

## 3. Canonical multi-seed results: RESOLVED

Main text and Supplementary Tables 4 and 5 report the multi-seed convention, and the
multi-seed tables are committed (`definitive_summary.csv`, `unified_multiseed.csv`).
`results/canonical/README.md` documents the dual convention and states which table backs
which number; the single-draw values (~0.46) are cited only where the paper says single-draw.
Every headline PDS, CI and per-gene value reproduces from a tracked file.

Two residuals: the harness code is not committed (B2), and `unified_summary5.csv` sits in the
public repo labelled "an alternative multi-seed run that disagrees with definitive_summary.csv;
kept for provenance, not cited". The label is honest but a reviewer who opens it will find two
disagreeing multi-seed runs; consider moving it to `results/deprecated/`.

## 4. JAK1 model-limited verdict: full chain present

| Link | File | Value |
|---|---|---|
| identification resolvable | `pairwise_resolvability.csv` | 15.4% nearest-sibling, 78.5% of pairs |
| identification classifier | `classifier_two_sample.csv` | ident AUROC 0.873 (perm 0.503) |
| oracle above chance | `oracle_ceiling.csv` | 0.792 (0.742-0.838) |
| allele residual reproducible | `residual_decomp_measurement.csv` | residual oracle 0.892 |
| current models at chance | `best_model_pds.csv` | best in-house 0.517 |
| published models comparable space | `definitive_summary.csv` | all scored in the full harness |

PerturbNet's reduced-space PDS is never compared against a full-space oracle: its 0.54 is
reported in the full harness, and the subspace deflation is a separate control. The
Discussion already reads "JAK1 therefore defines an evaluable but unresolved regime, not
evidence that allele-level prediction is intrinsically impossible", which is the calibrated
form. Add the n caveat (W9).

## 5. GATA1 boundary: consistent everywhere

Detection marginal (`D_self/D_null` 0.88), identification absent (0% nearest-sibling, 3.8% of
pairs, classifier 0.52), replicate ceiling modestly above chance (0.572, CI 0.549-0.595),
97.6% un-rankable at native depth and 100% at matched 50 cells per half. The SI decision
table row reads "marginal / marginal / not evaluable / measurement-limited / GATA1". No
passage claims GATA1 is fully unmeasurable, and none claims it carries a clear allele signal.

One nuance worth preserving: `residual_decomp_measurement.csv` gives GATA1
`frac_allele_residual` 0.706, i.e. most of its variance IS allele-specific, but that residual
is not reproducible (residual oracle 0.588). Do not use `frac_gene_shared` to support the
"direction is gene-shared" claim; the correct support is residual Pearson-delta near zero plus
the gene-mean tie.

## 6. Pilot claim: all six sub-conditions pass

`results/pilot_validation/README.md` documents the design: pilot = first 50 cells, eval =
next 2T cells, disjoint; pilot features from 25 to 50 cells; the label is split-half S > W on
the disjoint eval cells, i.e. detection-level.

* disjointness: PASS (README + `pilot_validate.py` docstring + Methods)
* detection-level: PASS (now stated in abstract, Results and Discussion)
* per-dataset not pooled: PASS. The README carries its own red line, that the pooled
  mechanistic AUROC 0.978 is inflated by cross-dataset base-rate separation and drops to 0.71
  after within-dataset rank normalisation, and that GATA1's 0.978 is a single-positive
  (1/183) statistic. The manuscript reports neither.
* no fake AUROC for single-class datasets: PASS. TP53 and KRAS are `evaluable=False` with
  blank AUROC in `rankability_predictor_honest.csv` and are reported as not evaluable.
* "prospective" means disjoint cells, not a new wet-lab cohort: PASS, stated in Methods.
* AUROC 0.85 / 0.94 / 0.99, mean 0.93: PASS (0.852 / 0.938 / 0.988, mean 0.9260).
  Permutation null 0.49: PASS. Holds at a 25-cell pilot: PASS (mean 0.92).

One residual: the per-dataset training-free values "0.84 to 0.98" exist only in the README
prose. `pilot_validate.py` writes only the pooled `mechanistic_auroc`. Minor
ARTIFACT_MISSING; fix by writing the per-dataset mechanistic AUROCs into the summary JSON.

## 7. Atlas duplication: keep the sentence, compress it

Audit basis for the ruling. The atlas numbers in the split-half subsection and in the pilot
subsection are the **same** native-depth quantity from the **same** files
(`canonical_numbers.json`, `unrankable_canonical.json`). The split-half subsection adds only
the 95% CIs and n; the pilot subsection additionally carries the native-versus-matched-50
contrast, which the split-half subsection does not. So the earlier occurrence is a preview,
not an independent support, and deleting the whole sentence would remove the only statement
that the floor is not specific to allele data.

Adopt the proposed compression verbatim:

> The same detection-level floor was also present in external gene-level Perturb-seq atlases,
> motivating the prospective design analysis below.

and move the CIs to the pilot subsection, where the numbers do their work. Delete the forward
reference to Fig. 6f from the split-half subsection; it is the only forward figure reference
in that subsection and the numbers it points at are restated there in full.

---

# PASS list

Verified exactly against a tracked file, no action needed.

**Benchmark composition.** `data/allele_perturb_bench.csv`: 472 rows minus the 2 wild-type
rows = 470 variants (TP53 98, KRAS 92, GATA1 254, JAK1 26); `n_cells` sum = **321,043**
including wild-type and control cells; median per-variant depth 929 / 1000 / 354 / 104. All
four numbers exact.

**Direction versus discrimination.** 19 non-null predictors span PDS 0.487 to 0.517
(mean 0.4989) with every bootstrap CI crossing 0.50 (`definitive_summary.csv`).
Pearson-delta spans 0.5545 to 0.6471 (`results_v4_exttheta.csv`, fig2c definition).
Permutation null mean 0.5000, range 0.4948 to 0.5038 over exactly 340 method-split-gene
combinations, 3.24% with permutation p < 0.05 (`permutation_null_pds.csv`).
Residual Pearson-delta -0.078 to 0.031 (`residual_decomp_models.csv`); residual oracle
0.480 / 0.478 / 0.588 / 0.892 (`residual_decomp_measurement.csv`).
DE gradient 0.6523 / 0.2548 / 0.1498 (fig2e definition).
Per-gene dissociation gap 0.3074 / 0.1877 / 0.1287 / -0.0404 (fig2f definition).

**Exclusions.** 15 of exactly 17 gene-by-split combinations; GATA1 low-depth PDS 0.1859;
leakage-conservative refit gives train-only max 0.5000 and mean delta -0.0115
(`leakage_control_summary.csv`); representation PDS 0.464 to 0.475 (`rep_grid_summary.csv`);
allele-blind coverage 4/470 = 0.85% (`allele_blind_demo.csv`); external models scGen 0.49
(0.46-0.52), scVIDR 0.50 (0.48-0.53), Biolord 0.49 (0.47-0.52), CellFlow 0.51 (0.48-0.54),
PerturbNet 0.54 (0.51-0.57), with Pearson-delta 0.18 / 0.13 / 0.13 / 0.12 / 0.13.

**Measurement floor.** `D_self/D_null` 0.965 / 1.004 / 0.878 / 0.210
(`canonical_numbers.json`); identification AUROC 0.493 / 0.498 / 0.520 / 0.873 against
permutation 0.50 and detection-versus-wild-type 0.533 / 0.498 / 0.508 / 0.965
(`classifier_two_sample.csv`); pairs resolvable 0 / 0 / 3.8 / 78.5%, nearest-sibling
distances 0.23 to 0.46 against replicate noise near 0.97 (`pairwise_resolvability.csv`);
five distance families tabulated (`metric_family_floor.csv`); un-rankable 100 / 100 / 97.6 /
10% for the allele genes and Replogle 55.29% (CI 52.84-57.51, n = 1832), Adamson 14.58%,
Norman 3.39%, with matched-50 values 54.89 / 38.54 / 11.16% (`unrankable_canonical.json`).

**Ceilings and benchmark validity.** Oracle 0.485 (0.445-0.525) / 0.500 (0.440-0.566) /
0.572 (0.549-0.595) / 0.792 (0.742-0.838) (`oracle_ceiling.csv`); all-cells oracle 0.486 and
0.492 and perfect-prediction-versus-noisy-truth 0.488 and 0.490
(`reviewer_controls/oracle_sensitivity.csv`); best in-house model 0.51 / 0.54 / 0.50 / 0.52
(`best_model_pds.csv`); above a ceiling of 0.65 both recovery probabilities span 0.921 to
1.000 (`controlled_recovery.csv`); nine-dataset resolution Adamson / Norman / Replogle / VCC
1.000, sci-Plex 0.772, JAK1 0.990, GATA1 0.512 shallow against 0.9455 native, KRAS 0.018,
TP53 0.014 (`benchmark_resolution/summary.csv`).

**Triage.** LODO AUROC mean 0.9574 (Replogle 0.9702, Norman 0.9337, Adamson 0.9303,
GATA1 0.953, JAK1 1.000), TP53 and KRAS `evaluable=False` with no AUROC reported
(`rankability_predictor_honest.csv`).

**Methods parameters.** 2,000 bootstrap resamples, percentile method, 95% level (matches
`N_BOOT = 2000`, `SEED = 0` in the fig2b/fig2c scripts); 50 random seeds per split-half
configuration with PCA-50 at native maximum depth; subsample depths 50, 100, 150, 300
(matches `split_half_power_curve.csv`); DE computed at `n_sub = 300`
(`alleleperturb/evaluation/de_metrics.py`); label-permutation and wild-type-versus-wild-type
controls both at 0.50 (`classifier_two_sample.csv`, `wt_control_diagnostic.csv`).

**Leakage disclosure.** The SI states 47% / 0% / 10% / 37% / 54% residue sharing for the
random, positional, mechanistic, low-depth and compatibility splits; `residue_leakage.csv`
gives 0.471 / 0.000 / 0.099 / 0.371 / 0.543. Exact.

**Supplementary tables.** Table 4 (published models) and Table 5 (per-gene backbone) match
their source tables cell by cell, including all four oracle CIs, the pairs-resolvable and
rankable-at-native columns, and the best-model column. The Table 4 PerturbNet footnote
(0.52 / 0.52 / 0.52 / P = 0.13) matches `perturbnet_subspace_test.json` exactly.

**Discussion.** The limitations paragraph already carries the four-gene and gene/assay
confound, the pilot's detection-level restriction, the oracle-is-an-empirical-reference
caveat and the pseudobulk caveat. The JAK1 sentence is already in the calibrated
"evaluable but unresolved" form.

---

# Repository and software fix list

Derived from the audit, to be folded into the packaging step.

1. Commit the canonical harness that produces the nine unbacked tables in `results/canonical/`,
   or narrow the Code availability sentence (B2).
2. Add the permutation test to `scripts/analysis/subspace_test.py`, make it write its JSON,
   and remove the hard-coded `/data/boom/NUS/VCCompass` path (B3).
3. Make `results/reviewer_controls/am_analysis.py` write `am_summary.csv`; remove its
   hard-coded path and document where `AM_4genes.tsv` comes from (B4).
4. Make `results/pilot_validation/pilot_validate.py` write the per-dataset mechanistic AUROCs
   and read from the committed `results/pilot_validation/*_pilot.csv` instead of
   `/data/boom/NUS/pilot_validation/` (item 6 residual).
5. Move `results/canonical/unified_summary5.csv` to `results/deprecated/` so the public repo
   does not ship two disagreeing multi-seed runs side by side (item 3 residual).
6. Package the four framework diagnostics with a CLI and `pyproject.toml` (B5).
7. Cosmetic: `manuscript/latex/figures/ED_fig1.pdf` renders as Supplementary Figure 1;
   rename or relabel. Same class of issue as printed Figure 3 being built from `fig4.pdf`
   after the Results reorder.
