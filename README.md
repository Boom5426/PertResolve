<div align="center">

# PertResolve

### Measurement resolution for fine-grained perturbation prediction

<p>
  <a href="https://github.com/Boom5426/PertResolve/actions/workflows/tests.yml"><img alt="Tests" src="https://img.shields.io/github/actions/workflow/status/Boom5426/PertResolve/tests.yml?label=tests&logo=github"></a>
  <a href="https://www.python.org/"><img alt="Python 3.10+" src="https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white"></a>
  <a href="https://github.com/Boom5426/PertResolve/blob/main/LICENSE"><img alt="MIT License" src="https://img.shields.io/badge/license-MIT-2ea44f"></a>
  <a href="https://huggingface.co/datasets/Boom5426/PertResolve_Bench"><img alt="Hugging Face dataset" src="https://img.shields.io/badge/data-Hugging%20Face-FFD21E?logo=huggingface&logoColor=black"></a>
  <a href="manuscript/PertResolve_manuscript.pdf"><img alt="Manuscript PDF" src="https://img.shields.io/badge/manuscript-PDF-B31B1B?logo=adobeacrobatreader&logoColor=white"></a>
</p>

<p>
  <strong>Measure the resolution of a perturbation measurement before using it to rank models.</strong>
</p>

</div>

<p align="center">
  <img src="assets/fig1_overview.png" alt="PertResolve overview: measurement questions, benchmark scope, allele representations, and evaluation protocol" width="920">
</p>

<p align="center"><em>Figure 1. PertResolve separates what a measurement can resolve from what a model predicts. <a href="assets/fig1.pdf">Open the vector PDF.</a></em></p>

> A benchmark cannot reliably distinguish two models on a biological distinction that its own measurement cannot reproduce. PertResolve reports the measurement resolution alongside the model-ranking evidence.

## ✨ What is PertResolve?

PertResolve is the software, benchmark metadata, analysis harness and final
paper artifacts accompanying **“Measurement resolution constrains fine-grained
perturbation prediction.”** It keeps four empirical questions separate:

| Question | Meaning |
| --- | --- |
| **Detection** | Can a perturbation be distinguished from the reference population? |
| **Identification** | Can the measured response be distinguished from its nearest competitor? |
| **Split-half reproducibility reference** | How similar are disjoint measurements of the same perturbation? |
| **Model-ranking resolution** | Which known differences between predictors can this measurement order? |

The central object is not a model score alone. It is the relationship between
signal, within-condition sampling noise, nearest-competitor separation and the
candidate pool being scored.

## 📌 At a glance

| Release component | Contents |
| --- | --- |
| **PertResolve-Bench** | 470 protein-coding variant conditions, plus 2 wild-type reference rows, across TP53, KRAS, GATA1 and JAK1 |
| **Allele-resolved arm** | 321,043 cells represented in the benchmark metadata table |
| **Broader panel** | 31 perturbation configurations spanning 14 public resources |
| **Evaluation layer** | Detection, nearest-competitor identification, split-half reproducibility reference, PDS, direction recovery and model-ranking resolution |
| **Public data release** | Small metadata and demo data in this repository; larger benchmark assets on Hugging Face |

## 🚀 Quick start

The following path is completely self-contained. It does not download a large
single-cell atlas, require an account, or depend on a machine-specific data
directory.

```bash
git clone https://github.com/Boom5426/PertResolve.git
cd PertResolve

python -m venv .venv
source .venv/bin/activate                 # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev,bench]"

# Repository and estimator tests
pytest

# A 160-profile, 939-feature resolution example
python examples/run_resolution_demo.py
```

The demo prints a deterministic `ResolutionReport` with independent detection,
identification, split-half reproducibility-reference and model-ranking-resolution
summaries. No single scalar is promoted to a pass/fail benchmark verdict. The demo file is
[`data/demo/pertresolve_resolution_demo.npz`](data/demo/pertresolve_resolution_demo.npz);
its provenance and limits are documented in
[`data/demo/README.md`](data/demo/README.md).

Typical output begins like this:

```text
matrix: 160 cells x 939 features
resolution report
  4 perturbations at 8 cells per group
  detection        75.0% of perturbations clear their split-half noise
  identification   50.0% are separable from their closest competitor
  split-half ref   0.938 empirical reproducibility reference
  model-ranking resolution  ...
```

## 🧭 Use the resolution API

The core package works directly with a NumPy matrix and one perturbation label
per cell. The matrix can contain expression features, a reduced embedding or
another fixed representation. The estimator performs disjoint cell splits so
the prediction-side and evaluation-side measurements do not reuse cells.

```python
import numpy as np
from pertresolve.resolution import resolution_report

with np.load("data/demo/pertresolve_resolution_demo.npz", allow_pickle=False) as demo:
    X = demo["X"]
    labels = demo["labels"].astype(str)

report = resolution_report(
    X,
    labels,
    control="control",
    depth=8,
    n_seeds=4,
    seed=7,
    n_boot=200,
)

print(report.summary())
```

For an AnnData file, install the optional I/O dependencies and use the CLI:

```bash
python -m pip install -e ".[io]"
pertresolve-resolution path/to/data.h5ad \
  --perturbation-key perturbation \
  --control non-targeting \
  --depth 50 \
  --out /tmp/pertresolve-resolution-demo
```

The CLI writes a JSON summary and a per-perturbation window table when the
window calculation is enabled. Output paths are explicit so a test run cannot
silently overwrite the committed `results/` tree.

## 🧬 Data and benchmark release

The repository keeps the small, inspectable metadata needed to understand the
benchmark in `data/`. The larger release bundle is hosted separately:

**[Boom5426/PertResolve_Bench on Hugging Face](https://huggingface.co/datasets/Boom5426/PertResolve_Bench)**

The benchmark has four allele-resolved arms:

- **TP53** and **KRAS**, derived from A549 Perturb-seq;
- **GATA1**, a base-editing perturbation screen;
- **JAK1**, an scSNV-seq perturbation screen.

After removing the two wild-type reference rows, the current benchmark table
contains 98 TP53, 92 KRAS, 254 GATA1 and 26 JAK1 variant conditions. The wider
panel adds 31 resource/configuration combinations so that resolution can be
compared across perturbation types, readouts and experimental designs.

To download the public benchmark bundle:

```bash
python -m pip install huggingface_hub
python - <<'PY'
from huggingface_hub import snapshot_download

snapshot_download(
    repo_id="Boom5426/PertResolve_Bench",
    repo_type="dataset",
    local_dir="data_external/PertResolve_Bench",
)
PY
```

The small tracked file [`data/pertresolve_bench.csv`](data/pertresolve_bench.csv)
is the repository-facing benchmark metadata table. The Hugging Face bundle also
contains the larger allele-response arrays, embeddings, barcode-level source
tables and the supplementary dataset workbook. Large matrices are deliberately
not copied into Git history.

## 🔁 Reproduce the paper

There are two reproducibility levels, with different input requirements.

### Level 1: repository-only checks and demo

These commands need only the repository and the Python dependencies:

```bash
python -m pip install -e ".[dev,bench]"
pytest
python examples/run_resolution_demo.py
```

They verify the estimator, repository contracts, split-half utilities and the
complete public API path on the bundled demo matrix.

### Level 2: paper-level analysis tables

The paper analyses use processed single-cell matrices and model predictions that
are too large to place in this Git repository. Point the analysis scripts to an
explicit workspace; never rely on a hidden machine-local default:

```bash
export PERTRESOLVE_DATA=/absolute/path/to/processed-pertresolve-workspace
export PERTRESOLVE_ATLAS_DIR=/absolute/path/to/perturbation-atlases  # when needed
```

The workspace used by the evaluation scripts contains the processed allele
arrays, benchmark metadata, model predictions and the released embedding files.
The analysis scripts fail with a readable input error when a required asset is
missing.

Representative paper-level commands write to scratch directories and leave the
committed results untouched:

```bash
# Split-half reproducibility reference
python scripts/analysis/oracle_ceiling.py \
  --base "$PERTRESOLVE_DATA" \
  --out /tmp/pertresolve-oracle

# Canonical multi-seed model comparison
python scripts/analysis/score_definitive.py \
  --base "$PERTRESOLVE_DATA" \
  --out /tmp/pertresolve-score

# Top-k split-half recovery
python scripts/analysis/split_half_topk.py \
  --base "$PERTRESOLVE_DATA" \
  --out /tmp/pertresolve-topk

# Re-run the broader panel for one AnnData resource
python scripts/analysis/resolution_panel_v2.py \
  --h5ad /path/to/resource.h5ad \
  --name resource_name \
  --perturbation-key perturbation \
  --control non-targeting \
  --out /tmp/pertresolve-panel

# Rebuild the Supplementary Data 1 workbook from committed tables
python scripts/make_dataset_table.py \
  --out /tmp/supplementary_data_1_datasets.xlsx
```

The analysis generators and their contracts are indexed in
[`scripts/analysis/README.md`](scripts/analysis/README.md). Committed derived
tables live in [`results/canonical/`](results/canonical/) and
[`results/resolution_panel_v2/`](results/resolution_panel_v2/). The final
manuscript PDF, Supplementary Information PDF and Supplementary Data 1 workbook
are already included under [`manuscript/`](manuscript/).

The final figure-generation source bundle is intentionally withheld from this
public staging pass while the manuscript text is being finalized. The static
Figure 1 preview above is included for orientation; the source bundle will be
added after the manuscript is accepted for release.

Detection, identification, split-half reproducibility and model-ranking resolution
are reported as separate axes. In particular, a split-half reference above any
operating point is not by itself evidence that a dataset can rank models; the
model-ranking result must be inspected directly.

## 🗂️ Repository map

```text
PertResolve/
├── pertresolve/                  # installable estimator and public API
├── data/                         # benchmark metadata and runnable demo data
├── examples/                     # end-to-end commands that work after install
├── results/                      # committed derived tables and provenance
├── scripts/analysis/             # paper analysis generators
├── manuscript/                   # exactly: 2 PDFs + 1 Supplementary Data workbook
├── tests/                        # estimator and repository-contract tests
├── assets/                       # static README visual assets
└── pyproject.toml                # package metadata and optional dependencies
```

See [`PROJECT_STRUCTURE.md`](PROJECT_STRUCTURE.md) for the longer directory
map and [`data/README.md`](data/README.md) for accession-level data notes.

## 📚 Citation

If you use the software or benchmark, please cite the paper listed in
[`CITATION.cff`](CITATION.cff). The code is released under the MIT License.

If PertResolve is useful for your work, please consider leaving a ⭐ on GitHub
so others can discover the project too.

```bibtex
@article{Li2026PertResolve,
  title   = {Measurement resolution constrains fine-grained perturbation prediction},
  author  = {Li, Bo and Zhang, Chengyang and Li, Mengran and Zhang, Bob and Wang, Lin and Tang, Zhenchao and Liu, Jun and Liu, Chengliang and Wei, Chen and Yi, Yuhao and Lv, Jiancheng and Zhang, Yang},
  year    = {2026},
  note    = {PertResolve manuscript}
}
```

<div align="center">

**Measure first. Compare second.** 🔬

</div>
