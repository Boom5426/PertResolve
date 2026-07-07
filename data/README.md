# Data

## Included

- `allele_perturb_bench.csv` — 472 variants with θ₆ features, split assignments,
  and per-variant cell counts.

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

After downloading, run the preprocessing script to generate z-scored
expression arrays:

```bash
python scripts/preprocess.py --raw_dir /path/to/downloads --out_dir data/
```

This produces:
- `joint_arrays.npz` (TP53 + KRAS, 939 shared HVGs)
- `gata1_arrays.npz` (2,477 HVGs)
- `jak1_arrays.npz` (2,000 HVGs)
