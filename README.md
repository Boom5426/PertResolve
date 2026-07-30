<div align="center">

<h1>🧬 AllelePerturb</h1>

<h3>Can models predict the transcriptional effect of <em>individual</em> protein-coding variants?</h3>

<p><b>A benchmark for allele-resolution single-cell perturbation prediction</b></p>

<p>
  <a href="LICENSE"><img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-green.svg"></a>
  <img alt="Python" src="https://img.shields.io/badge/python-3.11%2B-blue.svg">
  <img alt="Cells" src="https://img.shields.io/badge/single%20cells-321%2C043-orange.svg">
  <img alt="Variants" src="https://img.shields.io/badge/variants-470-9cf.svg">
  <img alt="Genes" src="https://img.shields.io/badge/genes-TP53%20%C2%B7%20KRAS%20%C2%B7%20GATA1%20%C2%B7%20JAK1-lightgrey.svg">
  <a href="https://github.com/Boom5426/AllelePerturb/stargazers"><img alt="Stars" src="https://img.shields.io/github/stars/Boom5426/AllelePerturb?style=social"></a>
</p>

<p>
  <a href="#-overview">Overview</a> ·
  <a href="#-key-results">Key Results</a> ·
  <a href="#-installation">Install</a> ·
  <a href="#-quick-start">Quick Start</a> ·
  <a href="#-data">Data</a> ·
  <a href="#-evaluation-protocol">Metrics</a> ·
  <a href="#-citation">Cite</a>
</p>

</div>

---

> [!IMPORTANT]
> **TL;DR** — Across 20 predictors, every method recovers the *direction* of a variant's transcriptional effect (Pearson-δ **0.60–0.68**), but **none reliably distinguishes individual alleles of the same gene** (PDS **0.49–0.52**, indistinguishable from chance 0.50). A split-half analysis traces this to a **measurement-resolution floor** set by effect size, sampling noise and sequencing depth — reframing variant-level benchmarking as a *power-aware* problem.

---

## ✨ Overview

Most single-cell perturbation benchmarks define a perturbation at the level of a **gene**, **drug** or **condition** ("knock down *TP53*", "apply compound X"). But many disease mechanisms are **allele-specific**: `TP53 R175H` unfolds the protein, `R273C` keeps it folded yet DNA-binding-dead, and `R248Q` shows dominant-negative / gain-of-function behavior.

**AllelePerturb** asks the finer question: *can a model tell these alleles apart from their single-cell transcriptional response?* It packages the publicly available single-cell perturbation data that resolve individual protein-coding variants into one benchmark, with a fixed evaluation protocol and generalization splits.

<div align="center">

| 🧫 Genes | 🔬 Variants | 🧮 Cells | 🧪 Technologies | 🎲 Splits | 📐 Metrics |
|:--:|:--:|:--:|:--:|:--:|:--:|
| 4 | 470 | 321,043 | 3 | 6 | 10 |

*TP53 · KRAS · GATA1 · JAK1   —   Perturb-seq · base editing · scSNV-seq*

</div>

---

## 🔑 Key Results

The central finding is a clean **direction–discrimination dissociation**:

<div align="center">

| | 🧭 **Direction recovery** | 🎯 **Allele discrimination** |
|:--|:--:|:--:|
| **Metric** | Pearson-δ | PDS (perturbation discrimination score) |
| **Result** | ✅ **0.60 – 0.68** | ❌ **0.49 – 0.52**  (chance = 0.50) |
| **Meaning** | models capture the shared, gene-level program | models cannot tell one allele from another |

</div>

> [!NOTE]
> **Why?** A split-half analysis shows the culprit is not the models but the **ground truth**: within-variant replicate noise approaches the variant-to-wild-type signal. When the *effect-size-to-noise window* is narrow, no method — however expressive — can rank held-out variants. The same floor appears in gene-level Perturb-seq atlases, where a large fraction of perturbations are un-rankable at native depth.

**Takeaway for practitioners:** before benchmarking models on variant-level data, first ask whether the ground truth is *measurable enough* to rank them. AllelePerturb provides the diagnostics and effect-size-conditioned guidance to do so.

---

## 📦 Installation

```bash
git clone https://github.com/Boom5426/AllelePerturb.git
cd AllelePerturb

conda create -n alleleperturb python=3.11 -y
conda activate alleleperturb
pip install -r requirements.txt
```

**Core dependencies:** `numpy` · `pandas` · `scipy` · `scikit-learn` · `matplotlib` · `seaborn` · `scanpy` · `fair-esm` · `torch`

---

## 🚀 Quick Start

Every figure is reproducible **from the pre-computed tables in `results/`** — no raw-data download required.

```bash
cd scripts/figures

python draw_fig2.py   # Direction–ranking dissociation
python draw_fig3.py   # Split-half measurement window
python draw_fig4.py   # Robustness across splits and metrics
python draw_fig5.py   # Rankability prediction + power-aware workflow
```

<details>
<summary><b>Re-run the full evaluation grid from raw data</b></summary>

```bash
# after downloading the raw arrays into data/ (see Data below)
python scripts/run_all_splits.py \
    --data_dir data/ \
    --output   results/results_v4_10metrics.csv
```

Runs 20 methods × 5 splits × 4 genes × 10 metrics (~30 min on a 24-core CPU).

</details>

---

## 🗂️ Repository Structure

```
AllelePerturb/
├── alleleperturb/          # Python package: bench loading, θ features, evaluation
├── data/                   # Benchmark tables (θ features, split assignments)
├── results/                # Pre-computed result tables (reproducibility backbone)
├── scripts/
│   ├── figures/            # Self-contained figure-drawing scripts
│   ├── baselines/          # Baseline method implementations
│   └── run_all_splits.py   # Full evaluation-grid runner
├── figures/                # Generated figures (composites + panels)
├── manuscript/             # LaTeX manuscript + references
├── requirements.txt
└── README.md
```

---

## 📊 Data

### Included in the repository

| File | Description |
|------|-------------|
| `data/allele_perturb_bench.csv` | 470 protein-coding variants (+2 WT reference rows): gene, protein, θ₆ biophysical features, 6 split assignments |
| `results/results_v4_exttheta.csv` | Full evaluation grid (canonical, external-θ) |
| `results/rankability_predictor_honest.csv` | Per-perturbation rankability predictor (leave-one-dataset-out) |

### Raw single-cell data (download separately, ~2 GB)

Needed only to **re-run the evaluation grid** — not to draw figures.

| Genes | Source | Assay |
|-------|--------|-------|
| **TP53 + KRAS** | [GEO GSE161824](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161824) | Perturb-seq (A549) |
| **GATA1** | [🤗 cyclopeta/PerturbNet_reproduce](https://huggingface.co/datasets/cyclopeta/PerturbNet_reproduce) | Base editing (HSPC) |
| **JAK1** | [Zenodo 10418435](https://doi.org/10.5281/zenodo.10418435) · ENA PRJEB48915 | scSNV-seq (HT-29) |

<details>
<summary><b>Download commands & preprocessing</b></summary>

```bash
cd data/

# TP53 + KRAS (Ursu et al. 2022)
for G in TP53 KRAS; do
  for F in processed.matrix.mtx.gz processed.genes.csv.gz variants2cell.csv.gz; do
    wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161824/suppl/GSE161824_A549_${G}.${F}
  done
done
python ../scripts/preprocess_ursu.py     # -> joint_arrays.npz

# GATA1 (Yu & Welch 2025)
python -c "from huggingface_hub import hf_hub_download; \
  hf_hub_download('cyclopeta/PerturbNet_reproduce','GATA1_standard_hvg_pert_filtered.h5ad',local_dir='.')"
python ../scripts/preprocess_gata1.py    # -> gata1_arrays.npz

# JAK1 (Cooper et al. 2024)
wget https://zenodo.org/records/10418435/files/scSNPseq_data.zip && unzip scSNPseq_data.zip -d jak1/
python ../scripts/preprocess_jak1.py     # -> jak1_arrays.npz  (needs R: scran + SingleCellExperiment)

# ESM-1v embeddings for all 470 variants
python ../scripts/extract_esm_embeddings.py --bench allele_perturb_bench.csv --output esm1v_embeddings.npz
```

</details>

---

## 📏 Evaluation Protocol

### Metrics

| Metric | Question | Range |
|--------|----------|:-----:|
| **PDS** *(primary)* | Is a prediction closer to its own target than to other variants? | 0–1 · **0.5 = chance** |
| **Pearson-δ** | Do predicted and true perturbation directions agree? | −1 to 1 |
| **DE overlap** | Fraction of top-50 DE genes shared | 0–1 |
| **Direction agreement** | Sign concordance across genes | 0–1 |
| **DE-LFC Spearman** | Rank correlation of log-fold-changes | −1 to 1 |
| **MAE** | Mean absolute error of predicted profiles | ≥ 0 |

*PDS is computed under cosine, L1 and L2 distances; **PDS-cosine** is the primary metric.*

### Generalization splits

| Split | Held out | Tests |
|-------|----------|-------|
| **Random** | 35% of variants | Standard generalization |
| **OOD-Position** | C-terminal half | Positional extrapolation |
| **OOD-Mechanism** | Hotspot / functional residues | Mechanistic extrapolation |
| **Cross-Gene** | one whole gene | Gene transfer |
| **Low-N** | variants < 200 cells | Low-depth regime |
| **Compatibility** | matched holdouts | Comparison with prior work |

---

## 📄 Manuscript

The LaTeX manuscript lives under `manuscript/latex/` and compiles with:

```bash
cd manuscript/latex && make
```

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

## 📝 License & Contact

Released under the [MIT License](LICENSE).
Questions and contributions welcome — please open an [issue](https://github.com/Boom5426/AllelePerturb/issues).

**Bo Li** · University of Florida, Department of Biomedical Engineering

<div align="center">
<sub>If AllelePerturb is useful for your work, consider leaving a ⭐ — it helps others find it.</sub>
</div>
