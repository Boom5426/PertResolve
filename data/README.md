# Data

## Included

- `pertresolve_bench.csv` — 472 rows: 470 protein-coding variant conditions plus 2 WT
  reference rows, with θ₆ features, split assignments and per-row cell counts.
- `demo/pertresolve_resolution_demo.npz` — a self-contained 160-profile,
  939-feature example for the public resolution API. It runs without external
  data; see [`data/demo/README.md`](demo/README.md).

## Full public release

The larger public benchmark bundle, including compact processed allele-response
products, embeddings, barcode-level source tables and the Supplementary Data
workbook, is available at
[Boom5426/PertResolve_Bench](https://huggingface.co/datasets/Boom5426/PertResolve_Bench).
The GitHub repository keeps the small metadata table and runnable demo locally so
the quick-start path remains lightweight. The Hugging Face README is the authority
for the exact files and their provenance.

## External data (download separately)

### TP53 + KRAS (Ursu et al. 2022, A549 Perturb-seq)

```bash
# From GEO GSE161824
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161824/suppl/GSE161824_A549_TP53.*
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161824/suppl/GSE161824_A549_KRAS.*
```

### GATA1 (PerturbNet, Yu & Welch 2025, HSPC base editing)

```bash
# From HuggingFace
pip install huggingface_hub
python -c "
from huggingface_hub import hf_hub_download
hf_hub_download('cyclopeta/PerturbNet_reproduce',
                'GATA1_standard_hvg_pert_filtered.h5ad',
                local_dir='.')
"
```

### JAK1 (Cooper et al. 2024, HT-29 scSNV-seq)

```bash
# From Zenodo
wget https://zenodo.org/records/10418435/files/scSNPseq_data.zip
unzip scSNPseq_data.zip
```

## Preprocessing

After downloading source data, run the preprocessing script to generate the
processed expression arrays used by the allele-resolved analyses:

```bash
python scripts/preprocess_gse161824_barcoded.py --help
```

The PertResolve Hugging Face release does redistribute the compact derived files
listed in its README, including ESM-1v embeddings, `real_deltas.npz`, the
processed JAK1 array and TP53/KRAS barcode and annotation tables. It does not
mirror the full cell-by-gene expression matrices or relicense the original GEO,
ENA, Zenodo or PerturbNet source data. Analysis scripts that need larger local
matrices accept them through `--base` or the `PERTRESOLVE_DATA` environment
variable; see the repository README.
