<div align="center">

<h1>🧬 AllelePerturb</h1>

<h3>Can your perturbation dataset arbitrate a model comparison at all?</h3>

<p><b>Measurement-resolution diagnostics for single-cell perturbation benchmarks, and the allele-resolution benchmark that motivated them</b></p>

<p>
  <a href="LICENSE"><img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-green.svg"></a>
  <img alt="Python" src="https://img.shields.io/badge/python-3.10%2B-blue.svg">
  <img alt="Dependencies" src="https://img.shields.io/badge/core%20deps-numpy%20%2B%20pandas-brightgreen.svg">
  <img alt="Tests" src="https://img.shields.io/badge/tests-23%20passing-success.svg">
  <img alt="Cells" src="https://img.shields.io/badge/single%20cells-321%2C043-orange.svg">
  <img alt="Variants" src="https://img.shields.io/badge/variants-470-9cf.svg">
  <a href="https://github.com/Boom5426/AllelePerturb/stargazers"><img alt="Stars" src="https://img.shields.io/github/stars/Boom5426/AllelePerturb?style=social"></a>
</p>

<p>
  <a href="#-run-it-on-your-data">Run it</a> ·
  <a href="#-what-governs-whether-a-benchmark-works">The finding</a> ·
  <a href="#-the-allele-benchmark">Benchmark</a> ·
  <a href="#-installation">Install</a> ·
  <a href="#-reproducing-the-paper">Reproduce</a> ·
  <a href="#-citation">Cite</a>
</p>

</div>

---

> [!IMPORTANT]
> **A leaderboard is only as good as the measurement under it.** If a dataset's ground truth
> cannot tell two perturbations apart, it cannot order two models that differ in how well
> they tell them apart, and its ranking reports the measurement as much as the methods. That
> is a property of the data, measurable before any model is trained. This repository makes it
> a one-line check, and reports what happens when you run it: across six published
> perturbation screens, **not one supports identification of individual perturbations** at 50
> cells per group.

---

## 🚀 Run it on your data

```bash
pip install -e ".[io]"
alleleperturb-resolution data.h5ad --perturbation-key perturbation --control non-targeting
```

```python
from alleleperturb.resolution import resolution_report

report = resolution_report(X, labels, control="non-targeting", depth=50)
print(report.summary())
```

```
verdict: detectable
  perturbations differ from the control but only 1% are separable from their closest
  competitor, so a discrimination score has nothing to reward
  414 perturbations at 50 cells per group, 1643 excluded
  detectable    51.2% of perturbations clear their replicate noise
  identifiable   1.4% are separable from their closest competitor
  ceiling      0.865 attainable by a second measurement
  ordering     orders predictors differing by 0.035 in quality (P(full order) = 0.96)
```

Three ordered questions, and the verdict names the first that fails:

| level | question | why it is not the previous one |
|---|---|---|
| **detectable** | does each perturbation differ from the control by more than its own replicate noise? | |
| **identifiable** | is each perturbation separable from its *closest competitor*? | several perturbations can each depart from the control while remaining indistinguishable from one another |
| **benchmarkable** | can the panel as a whole order predictors of graded, known quality? | a large screen can pool many individually unresolved perturbations into a benchmark that still orders predictors |

Where more cells would help, the report says so; where they would not, it says that too, and
the design calculator returns nothing rather than a number when the pilot separation's lower
bound does not clear zero. No finite depth recovers a separation that is not there, and a
confident cell count extrapolated from noise reads as a plan.

**The core diagnostics depend on numpy and pandas alone.** Everything they compute is a
distance, a mean or a rank. Requiring a scientific stack before someone can find out whether
their dataset can arbitrate a comparison would put the answer out of reach at the point it is
most useful.

---

## 🔍 What governs whether a benchmark works

A natural guess, and the one an earlier version of this analysis made, is that discrimination
depends on the ratio of the average separation between perturbations to the sampling noise.
**It does not.** Holding that ratio fixed at 0.59 and changing only the geometry of the
configuration, the attainable discrimination ranges from **0.999 to 0.632**.

The reason is that a discrimination score is a ranking question, which only a perturbation's
*nearest* competitor can spoil, while an average over pairs is dominated by far-apart pairs
the ranking never has to resolve. Over 1,452 controlled configurations:

| candidate axis | Spearman with attainable discrimination |
|---|---|
| mean squared separation over pairs | 0.767 |
| **median nearest-competitor separation** | **0.961** |

On four real datasets alone the two are indistinguishable (0.965 and 0.881), which is why
this took controlled configurations to settle rather than more data.

Two further results that change how such numbers should be read:

- **Ordering recovery is not monotone in resolution.** It peaks near an attainable
  discrimination of 0.705 and falls above it, because two predictors that both reach the
  ceiling become mutually indistinguishable. Only the ceiling can be inverted into a depth
  requirement.
- **Detection is common; identification is rare.** Across six screens run through the frozen
  criterion, detectable fractions span 26% to 90% while identifiable fractions span 0.5% to
  43%. The two largest gene-level atlases are 68.7% and 51.2% detectable against 0.5% and
  1.4% identifiable, and both still order graded predictors at probability 0.997 and 0.96.

Full derivations and the negative result that started them:
`docs/RESULT_COLLAPSE_REFUTED_2026-08-04.md`, `docs/RESULT_RESOLUTION_LAW_2026-08-04.md`,
`docs/RESULT_RESOLUTION_PANEL_2026-08-04.md`. The panel was
[pre-registered](docs/PREREG_RESOLUTION_PANEL_v1.md) before it was run, and the one
prediction that failed is reported as such.

```bash
pip install -e ".[dev]" && pytest      # 23 tests, simulations whose answer is known
```

The tests check the estimators against data with a known between-perturbation separation,
including that an unbiased estimate is recovered when the true separation is exactly zero.

---

## 🧪 The allele benchmark

The diagnostics came out of a concrete question. Most perturbation benchmarks define a
perturbation as a **gene**, **drug** or **condition**. Many disease mechanisms are
**allele-specific**: `TP53 R175H` unfolds the protein, `R273C` keeps it folded yet
DNA-binding-dead, `R248Q` behaves dominant-negatively. Can a model tell these apart from
single-cell transcriptional response?

AllelePerturb packages the publicly available single-cell data that resolve individual
protein-coding variants into one benchmark, with a fixed protocol: **470 protein-coding
variants, 321,043 cells**, ten metrics and six generalization splits.

The result is a clean **direction-discrimination dissociation**:

<div align="center">

| | 🧭 **Direction recovery** | 🎯 **Allele discrimination** |
|:--|:--:|:--:|
| **Metric** | Pearson-δ | PDS (perturbation discrimination score) |
| **Result** | ✅ **0.55 to 0.65** | ❌ **0.49 to 0.52** (chance = 0.50) |
| **Meaning** | models capture the shared, gene-level programme | models cannot tell one allele from another |

</div>

> [!NOTE]
> **The culprit is not the models.** A replicate ceiling, the discrimination a second
> measurement of the *same* variant attains, sits at chance for two of the four datasets: the
> measurement cannot reward even an essentially correct answer. For one dataset the ceiling
> is high while every model stays near chance, which is a genuine computational gap rather
> than a measurement one. Separating those two regimes is what the diagnostics above are for.

Scoring after removing the gene-shared programme does not rescue the models. That axis does
not have chance at 0.50 either: an uninformative prediction scores 0.524 on it, matching its
own permutation null exactly. Read against that null, six of 25 methods exceed it by 0.011 to
0.025 and none survives correction for the number tested, while the measurement's own ceiling
on the same axis rises, **widening** the gap between what the data permit and what the models
reach.

---

## 📦 Installation

```bash
git clone https://github.com/Boom5426/AllelePerturb.git
cd AllelePerturb
pip install -e .                 # resolution diagnostics; numpy and pandas only
pip install -e ".[io,bench]"     # add .h5ad reading and the benchmark analysis scripts
```

Reading `.h5ad` and reducing dimension need `anndata` and `scikit-learn` (`[io]`); the
benchmark loader and the analysis scripts under `scripts/` need `scipy` and `pyyaml`
(`[bench]`); `pytest` is `[dev]`.

---

## 🔁 Reproducing the paper

Two analyses read **only committed tables** and run on a fresh checkout with no external data:

```bash
python results/pilot_validation/pilot_validate.py --out /tmp/ap_out
python scripts/figures/rankability_predictor.py   --out /tmp/ap_out/rankability.csv
```

Both reproduce their committed counterparts in `results/` byte for byte.

Manuscript figures are built per panel under `manuscript/figures/figN/`, where each
`figN*_<panel>.py` writes a panel PDF and `figN_assemble.tex` composes them.
`manuscript/figures/check_panels.py` gates all 73 panels on size, type floor, overflow, row
width, per-figure height budget, single font family and label collisions.

### Locating data that is not in the repository

Analyses that touch single-cell data need two directories this repository does not ship,
because the files in them are either too large or not ours to redistribute. Each is named
explicitly, as an argument or an environment variable, and is never guessed:

| What | Argument | Environment variable | Must contain |
|---|---|---|---|
| **Processed allele data** | `--base` | `ALLELEPERTURB_DATA` | the per-gene cell matrices (`joint_arrays.npz` for TP53 and KRAS, `gata1_arrays.npz`, `jak1_arrays.npz`), the variant table `allele_perturb_bench.csv`, the protein embeddings `esm1v_embeddings.npz`, and the shared scorer at `unified/harness.py` |
| **Public perturbation atlases** | `--atlas-dir`, `--prep` | `ALLELEPERTURB_ATLAS_DIR` | the Replogle, Norman, Adamson, sci-Plex and Virtual Cell Challenge `.h5ad` files, under their original filenames |

Neither is inferred. When a location is missing, the script names both the argument and the
environment variable; when an input inside it is missing, the script prints the full path it
expected.

```bash
export ALLELEPERTURB_DATA=/path/to/processed-data
python scripts/run_all_splits.py --out /tmp/grid          # or: --base /path/to/processed-data
```

`--out` is required everywhere and is refused if it resolves inside `results/`, so a re-run
cannot overwrite a committed canonical table.

<details>
<summary><b>Script interfaces</b></summary>

| Script | Interface | Writes |
|---|---|---|
| `scripts/run_all_splits.py` | `[--base] --out` | `results_v4_10metrics.csv` |
| `scripts/analysis/oracle_ceiling.py` | `[--base] --out` | `oracle_ceiling.csv` |
| `scripts/analysis/pairwise_resolvability.py` | `[--base] --out` | `pairwise_resolvability.csv` |
| `scripts/analysis/classifier_two_sample.py` | `[--base] --out` | `classifier_two_sample.csv` |
| `scripts/analysis/metric_floor.py` | `[--base] --out` | `metric_family_floor.csv` |
| `scripts/analysis/controlled_predictors.py` | `[--base] --out` | `controlled_recovery.csv` |
| `scripts/analysis/floor_law_fit.py` | `[--base] --out` | `floor_law.csv` (superseded) |
| `scripts/analysis/resolution_scaling.py` | `[--base] --out` | `floor_law_v2.csv` |
| `scripts/analysis/resolution_sweep.py` | `[--base] --out` | `resolution_sweep.csv` |
| `scripts/analysis/resolution_law.py` | `--sweep --out` | calibration and depth prescriptions |
| `scripts/analysis/residual_axis.py` | `[--base] --out` | both scoring axes for every method |
| `scripts/analysis/subspace_test.py` | `[--base] --out` | subspace-test JSON |
| `scripts/figures/pipeline_v2.py` | `[--base] --out` | `all_metrics.csv`, `summary.json` |
| `scripts/figures/run_split_half_power_analysis.py` | `[--base] --out` | split-half power curve |
| `scripts/figures/rankability_audit.py` | `[--prep] --out <dataset>` | `<dataset>_rankability.csv` |
| `scripts/figures/rankability_predictor.py` | `[--table] --out` | per-dataset predictor results |
| `results/reviewer_controls/oracle_sensitivity.py` | `[--base] --out` | `oracle_sensitivity.csv` |
| `results/reviewer_controls/residual_decomp.py` | `[--base] --out` | residual-decomposition tables |
| `results/reviewer_controls/esm2_extract.py` | `[--base] --out [--device]` | four ESM2 feature `.npz` |
| `results/pilot_validation/pilot_features.py` | `[--base] [--atlas-dir] --out` | `<DS>_pilot.csv` |
| `results/pilot_validation/pilot_validate.py` | `--out` | `pilot_validation_summary.json` |
| `results/benchmark_resolution/benchmark_resolution.py` | `--out [--atlas-dir] <name>` | `<name>.json` |
| `results/benchmark_resolution/allele_resolution.py` | `[--base] --out` | `allele_<gene>.json` |

`results/reviewer_controls/esm2_extract.py` downloads the ESM2-650M weights (about 2.5 GB)
into the PyTorch hub cache on first use and expects a CUDA device.

Some tables are shipped as results rather than rebuilt here. No committed script regenerates
`results/results_v4_exttheta.csv` (read by the Fig. 2 and Fig. 4 panels),
`results/benchmark_resolution/summary.csv` (Fig. 5f) or `results/rankability_sensitivity.csv`
(Supplementary Fig. 1). The adapters used to run the published models (scGen, scVIDR,
Biolord, CellFlow, CPA, PerturbNet) and the shared evaluation harness are not part of this
repository.

</details>

---

## 📊 Data

### Included in the repository

| File | Description |
|------|-------------|
| `data/allele_perturb_bench.csv` | 470 protein-coding variants (plus 2 WT reference rows): gene, protein, θ₆ biophysical features, 6 split assignments |
| `results/canonical/` | Every table behind a manuscript number, each annotated with the panels it supports and each with a committed generator |
| `results/results_v4_exttheta.csv` | Full evaluation grid (canonical, external-θ) |

### Raw single-cell data (download separately, about 2 GB)

Needed only to re-run the evaluation grid, not to draw figures.

| Genes | Source | Assay |
|-------|--------|-------|
| **TP53 + KRAS** | [GEO GSE161824](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161824) | Perturb-seq (A549) |
| **GATA1** | [🤗 cyclopeta/PerturbNet_reproduce](https://huggingface.co/datasets/cyclopeta/PerturbNet_reproduce) | Base editing (HSPC) |
| **JAK1** | [Zenodo 10418435](https://doi.org/10.5281/zenodo.10418435) · ENA PRJEB48915 | scSNV-seq (HT-29) |

<details>
<summary><b>Download commands and preprocessing</b></summary>

```bash
cd data/

# TP53 + KRAS (Ursu et al. 2022)
for G in TP53 KRAS; do
  for F in processed.matrix.mtx.gz processed.genes.csv.gz variants2cell.csv.gz; do
    wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161824/suppl/GSE161824_A549_${G}.${F}
  done
done

# GATA1 (Yu & Welch 2025)
python -c "from huggingface_hub import hf_hub_download; \
  hf_hub_download('cyclopeta/PerturbNet_reproduce','GATA1_standard_hvg_pert_filtered.h5ad',local_dir='.')"

# JAK1 (Cooper et al. 2024)
wget https://zenodo.org/records/10418435/files/scSNPseq_data.zip && unzip scSNPseq_data.zip -d jak1/
```

**Preprocessing.** The workspace named by `ALLELEPERTURB_DATA` is expected to contain the
following arrays. Only the first has a committed producer:

| File | Produced by |
|---|---|
| `joint_arrays.npz` (TP53 + KRAS) | `scripts/figures/pipeline_v2.py`, which reads `raw/GSE161824_A549_*` under the workspace |
| `gata1_arrays.npz` | not in this repository; read if present |
| `jak1_arrays.npz` (needs R: scran + SingleCellExperiment) | not in this repository; read if present |
| `esm1v_embeddings.npz` | not in this repository |
| `allele_perturb_bench.csv` | not in this repository |

`raw/` is treated as immutable and is only ever read.

</details>

---

## 📏 Evaluation protocol

| Metric | Question | Range |
|--------|----------|:-----:|
| **PDS** *(primary)* | Is a prediction closer to its own target than to other variants? | 0 to 1 · **0.5 = chance** |
| **residual-PDS** | The same, after removing the gene-shared programme | **chance is 0.524, not 0.5**; read against a permutation null |
| **Pearson-δ** | Do predicted and true perturbation directions agree? | −1 to 1 |
| **DE overlap** | Fraction of top-50 DE genes shared | 0 to 1 |
| **Direction agreement** | Sign concordance across genes | 0 to 1 |
| **DE-LFC Spearman** | Rank correlation of log-fold-changes | −1 to 1 |
| **MAE** | Mean absolute error of predicted profiles | ≥ 0 |

PDS is computed under cosine, L1 and L2 distances; **PDS-cosine** is primary.

| Split | Held out | Tests | Scored |
|-------|----------|-------|:--:|
| **Random** | 35% of variants | standard generalization | ✅ |
| **OOD-Position** | C-terminal half | positional extrapolation | ✅ |
| **OOD-Mechanism** | hotspot / functional residues | mechanistic extrapolation | ✅ |
| **Low-N** | variants under 200 cells | low-depth regime | ✅ |
| **Compatibility** | matched holdouts | comparison with prior work | ✅ |
| **Cross-Gene** | one whole gene | gene transfer | ❌ defined but not scored |

Cross-gene is defined and released so others can use it, but is not scored here: the four
datasets do not share a gene space, so a model trained on one cannot emit a profile in
another's coordinates.

---

## 📖 Citation

```bibtex
@article{li2026alleleperturb,
  title   = {Measurement and model limits of allele-specific single-cell
             perturbation prediction},
  author  = {Li, Bo},
  year    = {2026},
  note    = {Manuscript in preparation}
}
```

---

## 📝 License and contact

Released under the [MIT License](LICENSE).
Questions and contributions welcome, please open an
[issue](https://github.com/Boom5426/AllelePerturb/issues).

**Bo Li** · Department of Artificial Intelligence, University of Macau

<div align="center">
<sub>If AllelePerturb is useful for your work, consider leaving a ⭐, it helps others find it.</sub>
</div>
