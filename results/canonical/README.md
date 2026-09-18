# Canonical results and reproduction index

This directory is the provenance map for the final manuscript PDFs in
[`manuscript/`](../../manuscript/). The final main-text numbering is Figures
1–6; the Supplementary Information contains Supplementary Figures 1–2 and
Supplementary Tables 1–6. The CSV/JSON files here are committed result
artifacts. The committed CSV/JSON artifacts are treated as frozen numerical outputs; this
index does not recompute, overwrite or hand-edit their numeric values.

## Reading the artifacts

Some committed panel artifacts preserve historical column names because their
values are frozen. In those files, `detectable_fraction`,
`identifiable_fraction` and `replicate_ceiling` should be read as the separately
reported detection, identification and split-half reproducibility-reference
axes. A historical `verdict` column is not the public resolution API and must
not be used to reconstruct a forced ladder. The current generator emits the
independent names `detection_fraction`, `identification_fraction`,
`split_half_reproducibility_reference` and `model_ranking_p_correct_order`.

The split-half quantity is an empirical reproducibility reference, not a hard
ceiling or bound. A value above 0.65 alone does not establish benchmarkability;
model-ranking resolution is a separate empirical output. Detection is measured
against the reference condition and identification against nearby candidates;
neither is treated here as a logical prerequisite for the other.

Status labels in the index mean:

- **confirmed**: the artifact and generator are explicitly tied by the script
  contract or a committed provenance check;
- **static**: a final PDF asset exists, but its rendering source is not part
  of the public software repository;
- **mapping unconfirmed**: no public generator-to-panel provenance has been
  confirmed. The mapping is intentionally not inferred from a suggestive filename.

All commands below write to scratch directories. Set the external input root
explicitly when a command needs the large processed arrays:

```bash
export PERTRESOLVE_DATA=/absolute/path/to/processed-pertresolve-workspace
```

The generators reject output paths inside the committed `results/` tree.

The final manuscript contains Figures 1–6 and Supplementary Figures 1–2.
This index records the public numerical artifacts and confirmed generator
mappings behind those figures; it does not infer panel provenance where that
mapping has not been established. The public repository also ships
`assets/fig1.pdf` as the static overview used by the README. Final journal
figure-layout sources are not included.

## Figure/panel → table → generator → input → command

| Figure/panel | Table or artifact | Generator | Required input | Reproduction command | Status |
| --- | --- | --- | --- | --- | --- |
| Fig. 1a–g | `assets/fig1.pdf` | Not public in this release | Final static Figure 1 PDF | No public rendering command; open `assets/fig1.pdf` | static |
| Fig. 2a | Final static `fig2.pdf` exists in local staging; no panel-specific table confirmed | — | — | — | static; data mapping mapping unconfirmed |
| Fig. 2b | `results/canonical/definitive_summary.csv` | `scripts/analysis/regenerate_canonical_pds.py` | Processed allele arrays, `preds5`, benchmark manifest | `python scripts/analysis/regenerate_canonical_pds.py --base "$PERTRESOLVE_DATA" --out /tmp/pertresolve-pds` | final manuscript panel; panel→table mapping mapping unconfirmed |
| Fig. 2c | `results/canonical/definitive_summary.csv` | `scripts/analysis/regenerate_canonical_pds.py` | Same as Fig. 2b | Same command as Fig. 2b | final manuscript panel; panel→table mapping mapping unconfirmed |
| Fig. 2d | `results/reviewer_controls/residual_decomp_measurement.csv` | `results/reviewer_controls/residual_decomp.py` | Processed allele arrays and canonical benchmark inputs | `python results/reviewer_controls/residual_decomp.py --base "$PERTRESOLVE_DATA" --out /tmp/pertresolve-residual-decomp` | final manuscript panel; table generator confirmed |
| Fig. 2e | `results/canonical/residual_axis_models.csv`, `residual_axis_per_gene.csv`, `residual_axis_ceiling.csv` | `scripts/analysis/residual_axis.py` | Processed allele arrays and prediction archive | `python scripts/analysis/residual_axis.py --base "$PERTRESOLVE_DATA" --out /tmp/pertresolve-residual` | final manuscript panel; table generator confirmed |
| Fig. 2f | Final static `fig2.pdf` exists in local staging; no panel-specific table confirmed | — | — | — | static; data mapping mapping unconfirmed |
| Fig. 2g | `results/canonical/released_model_table3_summary.csv` | `scripts/analysis/released_model_table3_summary.py` | Committed `unified_results5.csv` and `definitive_summary.csv` | `python scripts/analysis/released_model_table3_summary.py --check` | confirmed |
| Fig. 2h | `fig2_topk_summary.csv`, `fig2_topk_per_query.csv` | `scripts/analysis/fig2_topk.py` after `train_only_grid.py --dump-rank-counts` | Per-seed rank counts, canonical PDS table | `python scripts/analysis/fig2_topk.py --counts /tmp/pertresolve-grid/per_seed_rank_counts_trainonly.csv.gz --pds results/canonical/per_seed_variant_pds.csv.gz --out /tmp/pertresolve-topk` | final manuscript panel; table generator confirmed |
| Fig. 3a | Final static PDF exists in local staging; no numeric table confirmed | — | Local final figure asset only | — | static; public source pending |
| Fig. 3b | `oracle_ceiling.csv` and `definitive_summary.csv` | `oracle_ceiling.py` plus canonical model scorer | Processed arrays, predictions, benchmark manifest | `python scripts/analysis/oracle_ceiling.py --base "$PERTRESOLVE_DATA" --out /tmp/pertresolve-oracle`; then the Fig. 2b command | final manuscript panel; panel→table assembly mapping unconfirmed |
| Fig. 3c | Final static `fig3.pdf` exists in local staging; no panel-specific table confirmed | — | — | — | static; data mapping mapping unconfirmed |
| Fig. 3d | `results/canonical/pairwise_resolvability.csv` | `scripts/analysis/pairwise_resolvability.py` | Processed allele arrays and `pertresolve_bench.csv` | `python scripts/analysis/pairwise_resolvability.py --base "$PERTRESOLVE_DATA" --out /tmp/pertresolve-pairs` | final manuscript panel; table generator confirmed |
| Fig. 3e | `results/candidate_and_power/resolved_pair_accuracy.csv` | `scripts/analysis/resolved_pair_accuracy.py` | Processed allele arrays and prediction archive | `python scripts/analysis/resolved_pair_accuracy.py --base "$PERTRESOLVE_DATA" --out /tmp/pertresolve-pair-accuracy` | final manuscript panel; table generator confirmed |
| Fig. 3f | `results/canonical/split_half_topk.csv`, `split_half_topk_per_variant.csv` | `scripts/analysis/split_half_topk.py` | Processed allele arrays and canonical per-variant reference | `python scripts/analysis/split_half_topk.py --base "$PERTRESOLVE_DATA" --out /tmp/pertresolve-topk-reference` | final manuscript panel; table generator confirmed |
| Fig. 3g | Final static `fig3.pdf` exists in local staging; no panel-specific table confirmed | — | — | — | static; data mapping mapping unconfirmed |
| Fig. 4a | Final static `fig4.pdf` exists in local staging; no numeric table confirmed | — | Local final figure asset only | — | static; data mapping mapping unconfirmed |
| Fig. 4b | `results/canonical/fig4b_selfnull_dist.csv` | No single public panel generator confirmed | Native-depth window outputs | — | final manuscript panel; panel provenance mapping unconfirmed |
| Fig. 4c–d | `results/canonical/canonical_numbers.json`, `unrankable_canonical.json` | No single public panel generator confirmed | Native-depth window and native-depth eligibility summaries | — | final manuscript panel; data mapping mapping unconfirmed |
| Fig. 4e | `results/canonical/resolution_sweep.csv` | `scripts/analysis/resolution_sweep.py` | Processed allele arrays and controlled-geometry inputs | `python scripts/analysis/resolution_sweep.py --base "$PERTRESOLVE_DATA" --out /tmp/pertresolve-sweep` | final manuscript panel; table generator confirmed |
| Fig. 4f | `results/canonical/metric_family_floor.csv` | `scripts/analysis/metric_floor.py` | Processed allele arrays | `python scripts/analysis/metric_floor.py --base "$PERTRESOLVE_DATA" --out /tmp/pertresolve-metric-floor` | final manuscript panel; table generator confirmed |
| Fig. 5a–f | `results/canonical/resolution_panel_v2_table.csv` and `resolution_panel_v2_sensitivity.csv` | `scripts/analysis/build_resolution_panel_table.py` | Panel JSON runs plus v1 `resolution_panel.csv` | `python scripts/analysis/build_resolution_panel_table.py --panel results/resolution_panel_v2 --v1 results/canonical/resolution_panel.csv --out /tmp/pertresolve-panel-table` | final manuscript panel; aggregation confirmed |
| Fig. 6a | No panel-specific table confirmed; design schematic | — | — | — | mapping unconfirmed |
| Fig. 6b | `results/fig6_derived/prospective_auc.csv` | No public generator-to-panel mapping confirmed | Pilot/evaluation derived tables | — | final manuscript panel; data mapping mapping unconfirmed |
| Fig. 6c | `results/candidate_and_power/candidate_matched_summary.csv` | `scripts/analysis/candidate_matched_pds.py` | Processed allele arrays and prediction archive | `python scripts/analysis/candidate_matched_pds.py --base "$PERTRESOLVE_DATA" --out /tmp/pertresolve-candidate-pool` | final manuscript panel; table generator confirmed |
| Fig. 6d | `results/canonical/resolution_sweep.csv` | `scripts/analysis/resolution_sweep.py` | Controlled geometries and measured gene-by-depth inputs | See Fig. 4e command | final manuscript panel; table generator confirmed |
| Fig. 6e | `results/canonical/resolution_law.json`, `resolution_law_prescriptions.csv` | `scripts/analysis/resolution_law.py` | Scratch `resolution_sweep.csv` | `python scripts/analysis/resolution_law.py --sweep /tmp/pertresolve-sweep/resolution_sweep.csv --out /tmp/pertresolve-resolution-law` | final manuscript panel; table generator confirmed |
| Fig. 6f | `results/fig6_derived/depthmatch_summary.csv`, `depthmatch_paired.csv` | No public generator-to-panel mapping confirmed | Depth-matched atlas tables | — | final manuscript panel; data mapping mapping unconfirmed |
| Fig. 6g | `results/candidate_and_power/benchmark_power.csv` | `scripts/analysis/benchmark_power.py` | Processed allele arrays and controlled predictor inputs | `python scripts/analysis/benchmark_power.py --base "$PERTRESOLVE_DATA" --out /tmp/pertresolve-power` | final manuscript panel; table generator confirmed |
| Supplementary Fig. 1 | Final static `sfig1.pdf` exists in local staging; no panel-specific table confirmed | — | — | — | static; data mapping mapping unconfirmed |
| Supplementary Fig. 2 | Final static `sfig2.pdf` exists in local staging; no panel-specific table confirmed | — | — | — | static; data mapping mapping unconfirmed |
| Supplementary Table 1 | `data/pertresolve_bench.csv`, `data/hotspot_external_definition.txt` | No public table generator confirmed | Benchmark metadata and external annotation definition | — | mapping unconfirmed |
| Supplementary Table 2 | `configs/`, public API defaults and manuscript protocol | No public table generator confirmed | Configuration and fixed protocol definitions | — | mapping unconfirmed |
| Supplementary Table 3 | `results/canonical/released_model_table3_summary.csv` | `scripts/analysis/released_model_table3_summary.py` | Committed evaluation table and canonical PDS summary | `python scripts/analysis/released_model_table3_summary.py --check` | confirmed |
| Supplementary Table 4 | `results/canonical/oracle_ceiling.csv` plus model summaries | No single public table generator confirmed | Allele-resolved measurement and prediction summaries | — | mapping unconfirmed |
| Supplementary Table 5 | `results/canonical/resolution_panel_v2_table.csv` | `scripts/analysis/build_resolution_panel_table.py` | Panel JSON runs and v1 table | See Fig. 5 command | confirmed aggregation; table-to-PDF assembly pending |
| Supplementary Table 6 | `manuscript/supplementary_data_1_datasets.xlsx` and eligibility metadata | `scripts/make_dataset_table.py` for Supplementary Data 1 only | Committed metadata and result-file membership checks | `python scripts/make_dataset_table.py --out /tmp/supplementary_data_1_datasets.xlsx` | workbook generator confirmed; SI table mapping pending |

“Panel rendering pending” means the numeric artifact is traceable but the final
plotting source is not included. The public checkout intentionally does not
claim that the plot itself can be regenerated until that source bundle is
released.

## Primary artifact notes

- `definitive_summary.csv` is the primary multi-seed PDS summary. Its PDS and
  Pearson-related companion analyses use the paper's distinct held-out
  variant units and are not replaced by the generic resolution API.
- `oracle_ceiling.csv` is a historical filename for the split-half empirical
  reproducibility reference. It is not a perfect-score target and must not be
  interpreted as a universal upper limit.
- `controlled_recovery.csv` uses the generic controlled predictor family. Its
  `ceiling_pds` column is retained as a within-depth reference for that
  construction; it is not the paper's observed `ΔPDS` between model families.
- `resolution_panel_v2_table.csv` and the v1 panel artifacts preserve frozen
  legacy display fields. The independent public API and current aggregation
  generator are the semantic authority for new use.

The public repository is intended to reproduce the numerical analyses and
expose their provenance. Final journal figure-layout sources are not included;
the static Figure 1 asset used by the README is available at
[`assets/fig1.pdf`](../../assets/fig1.pdf).
