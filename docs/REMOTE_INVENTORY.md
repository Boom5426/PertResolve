# Remote Server Inventory (authoritative source of truth)

Compiled 2026-07-10 by systematic survey of the remote server. The local Ubuntu repo
(`/home/boom/ICLR/AllelePerturb`) is a copy; **the remote is authoritative** for all
compute, data, and model artifacts. This file maps what lives where, what is canonical,
and what is superseded/deprecated, so future sessions do not re-derive it.

Regenerate the raw manifest with:
`ssh 139.180.131.202 'find /data/boom/NUS/VCCompass /data/boom/VCData/VCCompass -type f -printf "%10s  %TY-%Tm-%Td  %p\n" | sort -k3'`

## 0. Connection and environments

- SSH: `ssh 139.180.131.202` (alias in `~/.ssh/config`: User `bob`, Port `2222`,
  key `~/.ssh/compute_key`). `/data/boom` -> symlink `/data/home_boom`.
- Conda env **`Agent`** (verified): python 3.11.14, numpy 2.4.6 / pandas 2.3.3 /
  scikit-learn 1.7.2 / scipy 1.14.1 / torch 2.10.0+cu126, RTX 4090 (24 GB). Used for the
  sklearn benchmark grid and general work.
- `unified/README.md` also references GPU envs **`pertbench`** (scGen/scVIDR/Biolord/PerturbNet),
  **`cellflow`** (CellFlow), **`cpa_legacy`** (CPA). NOTE: `conda env list` returned only
  base/Agent/Camo/NBE/SAM/boom, so these names must be reconciled before any re-run
  (they may be renamed, in another conda root, or venvs). The predictions already exist, so
  re-running is not required unless a method is changed.

## 1. Directory map (four project dirs)

| Dir | Role | Size | Git? |
|---|---|---|---|
| `/data/boom/NUS/VCCompass` | **active compute workspace**: run scripts, logs, the `unified/` harness, intermediate npz/checkpoints | 56 G, 609 files | not a git repo |
| `/data/boom/NUS/floor_audit` | **earlier per-method + cross-dataset floor workspace**: scGen/scVIDR/Biolord evals, GEARS demo, and the Replogle/Norman/Adamson rankability tables (section 6b) | 7.6 G | not a git repo |
| `/data/boom/NUS/methods` | **cloned external method repos + STATE runs**: CellFlow, STATE, scDFM source; STATE per-gene training (section 6c) | 12 G | contains vendored repos |
| `/data/boom/VCData/VCCompass` | **data + artifact store**: immutable `raw/`, big arrays, `.pt` checkpoints, mirrored `results/`, `summary.json` | 4.1 G, 53 files | not a git repo |

Shared files (`allele_perturb_bench.csv`, `summary.json`, arrays, checkpoints, `results/`)
are byte-identical across NUS/VCCompass and VCData/VCCompass (VCData is the curated mirror;
NUS is where work happens). Do not confuse with the unrelated `/data/boom/Agent/MVCBench` (a
different paper).

## 2. Immutable raw data (never modify, never write into `raw/`)

`raw/GSE161824_A549_{KRAS,TP53}.processed.matrix.mtx.gz` (816 M / 829 M, dated 2022-02-16),
plus `.genes.csv.gz` and `.variants2cell.csv.gz`. This is the Ursu 2022 A549 Perturb-seq
coding-variant screen (TP53 + KRAS). GATA1 and JAK1 come from external sources (see below),
not from this `raw/`.

## 3. Benchmark data + model artifacts

**Grid arrays (canonical model inputs, in both dirs):**
- `joint_arrays.npz` (386 M) = TP53 + KRAS (keys `Xtp/vtp/THtp`, `Xkr/vkr/THkr`).
- `gata1_arrays.npz` (658 M) = GATA1 base-editing screen (from PerturbNet's processed matrix
  + predefined holdout; HF repo `cyclopeta/PerturbNet_reproduce`). Keys `X`, `cell_variants`.
- `jak1_arrays.npz` (19 M) = JAK1 scSNV-seq, HT-29 + IFN-gamma (ENA PRJEB48915, Zenodo 10418435).
  Keys `X`, `variant_labels`. **JAK1 is the positive control** (D_self/D_null = 0.21).
- `esm1v_embeddings.npz` (2.3 M) = ESM-1v per-variant embeddings, keys `{GENE}__{variant}` incl `__WT`.
- `theta_v2.csv`, `gata1_theta.csv`, `jak1_theta.csv`, `wt_seqs.json` = biophysical theta features.

**"g1" flagship counterfactual arrays (GATA1, mostly superseded by unified):**
`g1_real_cells.npz` (511 M), `g1_cf_cells_flagship.npz` (910 M), `g1_deltas_{true,cf,recon}.npz`,
`g1_metadata.parquet`. WARNING per unified/README: `g1_real_cells.npz` GATA1 WT is anti-correlated
(cos -0.45) with the grid arrays; the unified harness deliberately does NOT use it.

**Model checkpoints (`.pt`, in-house models, in both dirs):**
`model_biophys.pt`, `model_dosage.pt`, `model_full.pt`, `model_random.pt` (~20 M each),
`model_cfm{,_v2,_v3}.pt` (~29 M, conditional flow-matching decoder iterations).

## 4. Bench / theta construction (provenance chain)

- `prepare_data_remote.py`, `fix_process.py` = raw mtx -> per-gene arrays.
- `g5_esm1v.py` = ESM-1v extraction (apply_mut). ESM keys must be `{GENE}__{variant}`.
- `theta_v2.csv` + per-gene theta CSVs = biophysical features (d_hydro/d_vol/d_charge/fold_core/
  cat_switch/is_hotspot) + `split1..6_role`. NOTE (from memory): the canonical theta z-score +
  split-role generator that emits `allele_perturb_bench.csv` is NOT fully archived on-server;
  rebuilding the 472-row bench bit-for-bit is a known prerequisite before adding genes.
- `allele_perturb_bench.csv` (472 variants) and `allele_perturb_bench_v2.csv` (adds
  `is_hotspot_leaked_OLD`) = the master bench tables. `all_splits_v4.py` reads v1; the exttheta /
  permutation runners read v2.

## 5. Main benchmark grid (legacy lineage, SUPERSEDED by `unified/` for method comparison)

- `all_splits_v4.py` (canonical legacy grid runner) + `all_splits_v4_exttheta.py` (de-leaked theta,
  external hotspot only). v1/v2/v3 (`all_splits_runner.py`, `all_splits_v2.py`, `all_splits_v3.py`)
  are earlier iterations.
- Outputs: `results_v4_10metrics.csv` (~10k rows, the paper's per-gene metric grid),
  `results_v4_exttheta.csv`, `results_all_splits_all_methods.csv`, `all_metrics.csv`, `summary.json`
  (the 4-baseline overall PDS ~0.70 table = the null-result headline).
- Legacy deep-model runners (SUPERSEDED by unified/run_*5.py): `cellflow_run{,_v2,_v3}.py`,
  `cpa_run{,_v2,_v3}.py`, `cfm_decoder{,_v2,_v3}.py` -> `cfm_eval{,_v2,_v3}.csv`.

## 6. `unified/` — CANONICAL method comparison (single source of truth, 2026-07-09/10)

One deterministic harness scoring **every** method on identical footing. THIS is what the
manuscript's main comparison should trace to, not the scattered legacy eval CSVs.

- `harness.py` = single source of truth: data loading, deterministic pseudobulk `real_deltas.npz`,
  tie-aware `PDS_cos` copied verbatim from `all_splits_v4.score_variant`, splits 1/2/3/5/6
  (split4 = cross-gene excluded, matching the paper), theta map. Also builds `grid_cbv.npz` (515 M,
  cached per-variant cells, 300-cap deterministic).
- Runners: `run_grid5.py` (18 sklearn heads {Ridge,Lasso,RF,GBoost,KNN,MLP}x{theta,ESM,ESM+theta}
  + Gene-mean + WT-null; Agent env), `run_scgen5.py` (scGen + scVIDR), `run_biolord5.py`,
  `run_perturbnet5.py` (cINN, ESM-conditioned), `run_cellflow5.py`, `run_cpa5.py` (did not converge).
- `preds5/*.npz` = per-method per-variant predicted deltas (auditable handoff format). Present:
  Biolord, CellFlow, PerturbNet, scGen, scVIDR + full 18-head grid + Gene-mean + WT-null.
  `variant-CPA.STATUS.txt` documents CPA non-convergence.
- Scoring: `score5_fast.py` -> `unified_summary5.csv` (+ per-variant `unified_results5.csv`);
  `score_definitive.py` (2026-07-10, newest) -> `definitive_summary.csv` (bootstrap CI);
  `score5_multiseed.py` -> `unified_multiseed.csv`. `subspace_test.py` = PerturbNet PCA-subspace
  artifact test. `allele_blind_demo.{py,json,csv}` = the GEARS/STATE category-error demonstration.

**Result (24/25 methods at chance, CI crosses 0.5):**
- Published SOTA run to completion, all at chance: scGen ~0.49, scVIDR ~0.50, Biolord ~0.52,
  CellFlow ~0.51. The 18 sklearn heads: 0.47-0.52.
- **PerturbNet ~0.54 is a PCA-subspace artifact** (`subspace_test.py`: a random vector confined to
  its 50-dim WT-PCA subspace already scores 0.516; effective subspace-chance ~0.52; empirical
  permutation null 0.538; direction recovery ~0, pearson-delta 0.13). Report with this caveat.
- **variant-CPA did not converge** (adversarial autoencoder -> NaN on continuous 6-dim theta;
  OOD embedding path is correct). Documented numerical-stability limitation.
- **STATE / GEARS / scGPT are allele-blind by construction** (gene token / GO node / one-hot
  perturbation identity): they cannot ingest a held-out variant's features, collapse all variants
  of a gene to one prediction (coverage 4/472 = 0.85%), so an allele-level PDS is a category error.
  Belong in a SEPARATE panel, NOT the PDS ranking. Empirical anchors: Gene-mean reference = 0.468
  (chance); STATE native discrimination_score on JAK1 ~0.556 independently reproduces the floor.

**Decided final table structure (unified/README.md, 2026-07):**
- Main: 4 representative baselines (Ridge-theta / GBoost-theta / MLP-theta / KNN-theta) + Gene-mean
  + WT-null references + variant-conditionable SOTA (scGen, scVIDR, Biolord, CellFlow), all at chance;
  PerturbNet with the subspace-artifact caveat.
- Supplementary: full 18x3 head grid (`unified_summary5.csv`).
- Separate "allele-blind by construction" panel: scGPT, GEARS, STATE (not in PDS ranking).
- "Did not converge": variant-CPA.

## 6b. `floor_audit/` — earlier per-method evals + cross-dataset floor (KEY for generalization)

`/data/boom/NUS/floor_audit/` (7.6 G, 2026-07-06) is the round that fed the unified harness and
the manuscript's generalization claim. Big scGen/scVIDR/Biolord training logs live under
`lightning_logs/` and `scvi_log/`. Result CSVs in `floor_audit/results/`:
- Per-method allele evals (earlier round; superseded for the main table by `unified/`, but the
  source of the numbers catalogued in memory [[benchmark-methods-audit]]): `biolord_eval.csv`,
  `groupB_scgen_eval.csv` (scGen + scVIDR), `groupA_matrix_eval.csv` (sklearn grid),
  `gears_allele_blindness.csv` + `gears_demo_summary.json` (GEARS allele-blindness demonstration).
- **Cross-dataset rankability / floor (`results/floor_audit/`) — the abstract's "same floor in
  gene-level atlases" evidence, and the asset for the Replogle generalization (review direction 3):**
  `Replogle_rankability.csv` (6.5 M) + `Replogle_floor_table.csv`, `Norman_rankability.csv` +
  `Norman_floor_table.csv`, `Adamson_rankability.csv` + `Adamson_floor_table.csv`; the four
  AllelePerturb genes `TP53/KRAS/GATA1/JAK1_rankability.csv`; plus `floor_vs_n_curves.csv`,
  `dataset_inventory.csv`, `s_vs_w_summary.csv`.

## 6c. `methods/` — vendored external method repos + STATE runs

`/data/boom/NUS/methods/` (12 G): cloned source for the external SOTA:
- `CellFlow/` (Theis lab OT/flow), `state/` (Arc Institute STATE; has MODEL_LICENSE / AUP),
  `scDFM/` (single-cell diffusion model).
- `state_data/` (per-gene `{KRAS,JAK1,TP53,GATA1}_dir` + `splits.json`) and `state_runs/`
  (`{GATA1,KRAS,TP53,JAK1,JAK1_smoke}/checkpoints/`, `eval_final.ckpt`): STATE was trained per gene
  here. STATE conditions on one-hot / gene-token perturbation identity, so it is allele-blind for
  truly held-out variants; its native discrimination_score on JAK1 (~0.556) reproduces the floor.
  STATE/GEARS therefore belong in the "allele-blind by construction" panel (section 6), NOT the PDS
  ranking. Note: these were RUN (contradicting any "never attempted" reading); the design decision
  is to exclude them from the allele-level ranking as a category error, with a quantitative anchor.

## 7. Debt / rigor analyses (leakage + pooling guards) -> `results/`

- `debt1_no_leakage_retrain.py` -> `results/debt1_no_leakage.csv`, `debt1_biophys_per_gene.csv`,
  `debt1_summary.csv`, `debt1_decision.json`: features masked at training time; VC_biophys PDS
  0.503 [0.470,0.537], CI crosses 0.5.
- `debt2_pds_statistics.py` -> `results/debt2_{pds_stats,permutation_null,pervariant_pds,pool_inflation,wilcoxon}.csv`,
  `debt2_decision.json`: permutation null + paired Wilcoxon + pool-inflation guard (perm-p ~0.09 n.s.).

## 8. Power / rankability / detection-limit -> `results/` and `results/power/`

- `run_split_half_power_analysis.py` -> `results/power/split_half_power_{curve,summary_by_gene,table}.csv`,
  `decision_after_split_half_floor.md`: the D_self/D_null measurement-floor analysis.
- `r3_detection_limit_power_curve.py` -> `results/r3_{critical_n,power_curve,required_n_prescription,
  underpowered_genes}.csv`, `r3_decision.json`: effect-size-conditioned required-N prescription.

## 9. PerturbNet-specific

- `external_PerturbNet/` = full shallow git clone of the PerturbNet repo (source + tutorials +
  GATA1 example notebooks). `perturbnet_canonical.py` + `perturbnet_canonical_stats.py` ->
  `results/perturbnet_canonical.{csv,preds.npz}`, `perturbnet_canonical_{decision,stats}.json`
  (heldout PDS 0.533, full-pool 0.572, pearson-delta 0.274).
- `deprecated_leaky/perturbnet_comparison.csv` (+ README) = **DO NOT USE**: the old `g5_cinn.py`
  PDS_within_gene_cosine 0.61 was a non-canonical, PCA-50, non-tie-aware, leakage-prone metric.
  Superseded by the canonical analyses above.

## 10. RUNX1 raw (blocked dataset-expansion attempt, 2026-07-09)

`runx1/raw/SRR29286119_{1,2}.fastq.gz`, `SRR29286122_{1,2}.fastq.gz` (~53 GB total raw fastq) +
`ena_manifest.tsv`, `fastq_urls.txt`, `parent_manifest.tsv`, `aria2.log`. RUNX1 (SEUSS) needs raw
SRA reprocessing since per-cell variant assignments are not deposited reusably. Integration is
paused (per memory, awaiting author-provided processed genotypes). These fastqs are large scratch.

## 11. Manuscript compile dirs on server (older copies; local is ahead)

`latex/` and `latex_full/` hold Jul-8 pre-edit `.tex` + `references.bib` + `figures/fig1-5.pdf`,
`ED_fig1.pdf`, compiled PDFs. Local `manuscript/latex/` is AHEAD (Data availability section, SI).
Direction = push local -> remote when a compute-side compile with canonical figures is needed;
back up the remote copy first.

## 12. What to trust (canonical vs superseded), quick reference

- Method comparison: **`unified/` (definitive_summary.csv + unified_summary5.csv)**. NOT the legacy
  `cfm_eval*.csv`, NOT `deprecated_leaky/`, NOT `perturbnet_comparison.csv`.
- Null-result headline (4 baselines ~0.70): `summary.json` + `results_v4_10metrics.csv`.
- Leakage/pooling guards: `results/debt1_*`, `results/debt2_*`.
- Power / required-N: `results/power/*`, `results/r3_*`.
- Raw data: `raw/` (immutable).

## 13. Known gaps / risks (verify before relying)

- Env names `pertbench`/`cellflow`/`cpa_legacy` not found by `conda env list`; reconcile before re-run.
- Canonical theta+split-role generator for `allele_perturb_bench.csv` not fully archived on-server.
- `definitive_summary.csv` (score_definitive, 7-10) vs `unified_summary5.csv` (score5_fast) give
  slightly different numbers (e.g. PerturbNet 0.542 vs 0.561; scVIDR NaN handling); decide one
  canonical scorer before quoting numbers in the manuscript.
- Fig 5 LODO rankability predictor (AUROC ~0.974, per memory) reproducibility not yet confirmed
  against a committed script.
