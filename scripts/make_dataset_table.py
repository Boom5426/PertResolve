#!/usr/bin/env python3
"""Build Supplementary Data 1: the dataset inventory for the PertResolve manuscript.

Scope. This workbook DESCRIBES the datasets: what they are, how large they are,
where they came from and which analysis arm each one enters. It deliberately
carries no per-dataset results. The measured quantities live in the
Supplementary Tables of the manuscript (Table 5 per-gene backbone and
reproducibility ceiling, Table 6 signal and noise, Table 7 depth prescriptions,
Table 9 the pre-registered panel), and duplicating them here would create a
second, drift-prone copy of numbers that are already typeset.

Numbering. Nature Portfolio numbers Supplementary Tables and Supplementary Data
as separate sequential series. The Supplementary Information PDF of this
manuscript already contains Supplementary Tables 1-10, so this standalone
workbook is Supplementary Data 1, not Supplementary Table 1.

Provenance. Two kinds of field, kept apart on purpose:
  * quantities measured in this study, read from the committed result files and
    the benchmark table, never typed in;
  * published metadata for the seven external screens, each verified against a
    primary source (publisher full text, repository record or bibliographic
    API) and carried with its URL on the Sources sheet.
A field that could not be verified from a source that was actually retrieved is
written as NOT_VERIFIED rather than filled in, so the workbook cannot be read as
asserting something this study did not check.

Sheets
------
Datasets            eleven datasets, description and scale
Analysis_roles      which dataset enters which analysis arm (the arms are not nested)
Column_definitions  what each column means, and where the results live instead
Sources             every external metadata value with its source URL, plus the
                    internal result files behind the measured columns

Usage
-----
    python scripts/make_dataset_table.py --out manuscript/supplementary_data_1_datasets.xlsx

``--out`` is required and has no default, matching the write guard the analysis
scripts already use: nothing here may silently overwrite a committed file.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
sys.path.insert(0, str(REPO))
from pertresolve.paths import reject_repo_results  # noqa: E402
UNVERIFIED = "NOT_VERIFIED"

# ---------------------------------------------------------------------------
# Primary allele-resolved screens. Assay, context and accession are as stated in
# the manuscript Methods (Benchmark construction, Data availability); the
# bibliographic fields are from references.bib.
# ---------------------------------------------------------------------------
ALLELE_STATIC = {
    "TP53": dict(
        protein_len_aa=393, uniprot="P04637",
        modality="Protein-coding variant (allele)", assay="Perturb-seq",
        cell_context="A549 (lung adenocarcinoma)", species="Homo sapiens",
        deposit="GEO GSE161824", obtained_from="GEO GSE161824",
        reference="Ursu et al. 2022, Nature Biotechnology 40:896-905",
        doi="10.1038/s41587-021-01160-7"),
    "KRAS": dict(
        protein_len_aa=189, uniprot="P01116",
        modality="Protein-coding variant (allele)", assay="Perturb-seq",
        cell_context="A549 (lung adenocarcinoma)", species="Homo sapiens",
        deposit="GEO GSE161824", obtained_from="GEO GSE161824",
        reference="Ursu et al. 2022, Nature Biotechnology 40:896-905",
        doi="10.1038/s41587-021-01160-7"),
    "GATA1": dict(
        protein_len_aa=413, uniprot="P15976",
        modality="Protein-coding variant (allele)", assay="Base editing (Perturb-BE-seq)",
        cell_context="Human haematopoietic stem and progenitor cells", species="Homo sapiens",
        deposit="GEO GSE215253",
        obtained_from="Processed matrix and holdout splits distributed with PerturbNet "
                      "(Hugging Face cyclopeta/PerturbNet_reproduce)",
        reference="Martin-Rufino et al. 2023, Cell 186:2456-2474.e24",
        doi="10.1016/j.cell.2023.03.035"),
    "JAK1": dict(
        protein_len_aa=1154, uniprot="P23458",
        modality="Protein-coding variant (allele)", assay="scSNV-seq",
        cell_context="HT-29 (colorectal), interferon-gamma stimulated", species="Homo sapiens",
        deposit="ENA PRJEB48915", obtained_from="Processed objects from Zenodo record 10418435",
        reference="Cooper et al. 2024, Genome Biology 25:20",
        doi="10.1186/s13059-024-03169-y"),
}

# ---------------------------------------------------------------------------
# The seven broader screens. Every value below was verified against the source
# recorded in EXTERNAL_SOURCES; nothing here is from recall. Where a paper
# reports no study-wide total, that is stated rather than summed by hand.
# ---------------------------------------------------------------------------
EXTERNAL_STATIC = {
    "Adamson": dict(
        modality="CRISPRi gene knockdown (dCas9-KRAB)",
        assay="Perturb-seq (single-cell CRISPR screen)",
        cell_context="K562 (chronic myeloid leukaemia)", species="Homo sapiens",
        pub_perturbations="91 sgRNAs targeting 82 genes (unfolded protein response screen)",
        pub_cells="Reported per experiment; no study-wide total stated",
        deposit="GEO GSE90546", obtained_from="scPerturb data portal",
        reference="Adamson et al. 2016, Cell 167:1867-1882.e21",
        doi="10.1016/j.cell.2016.11.048"),
    "Norman": dict(
        modality="CRISPRa, single and paired (SunTag)",
        assay="Perturb-seq (single-cell CRISPR screen)",
        cell_context="K562 (chronic myeloid leukaemia)", species="Homo sapiens",
        pub_perturbations="287 perturbations",
        pub_cells="~110,000 cells (median 273 per condition) in the pooled experiment",
        deposit="GEO GSE133344", obtained_from="scPerturb data portal",
        reference="Norman et al. 2019, Science 365:786-793",
        doi="10.1126/science.aax4438"),
    "Replogle": dict(
        modality="CRISPRi gene knockdown, genome-scale",
        assay="Perturb-seq (single-cell CRISPR screen)",
        cell_context="K562 and RPE1", species="Homo sapiens",
        pub_perturbations="9,866 expressed genes (K562 day-8 genome-scale screen)",
        pub_cells=">2.5 million cells",
        # The manuscript's Methods say "GEO" nowhere for this dataset, and the
        # paper's own Data availability names only SRA. Recording a GSE number
        # here would be wrong; see the Sources sheet.
        deposit="SRA BioProject PRJNA831566 (no GEO accession); processed data at "
                "figshare+ 10.25452/figshare.plus.20029387",
        obtained_from="scPerturb data portal",
        reference="Replogle et al. 2022, Cell 185:2559-2575.e28",
        doi="10.1016/j.cell.2022.05.013"),
    "VCC": dict(
        modality="CRISPRi gene knockdown (dual-guide dCas9-KRAB)",
        assay="Single-cell CRISPR screen (Virtual Cell Challenge training split)",
        cell_context="H1 human embryonic stem cells", species="Homo sapiens",
        pub_perturbations="150 target genes in the training split (of 300 total)",
        pub_cells="183,097 cells over the 150 training targets (sum of the official "
                  "pert_counts_Training.csv; excludes non-targeting controls)",
        deposit="Arc Virtual Cell Atlas, Google Cloud Storage bucket; no INSDC accession",
        obtained_from="Arc Virtual Cell Atlas",
        reference="Roohani et al. 2025, Cell 188:3370-3374",
        doi="10.1016/j.cell.2025.06.008"),
    "sciPlex": dict(
        modality="Small molecule, 188 compounds at four doses",
        assay="sci-Plex (sci-RNA-seq with nuclear hashing)",
        cell_context="A549, K562 and MCF7 in the main screen", species="Homo sapiens",
        pub_perturbations="188 compounds, each at four doses",
        pub_cells="649,340 cells",
        deposit="GEO GSE139944", obtained_from="GEO GSE139944",
        reference="Srivatsan et al. 2020, Science 367:45-51",
        doi="10.1126/science.aax6234"),
    "McFarland": dict(
        modality="Small molecule (13 drugs) and CRISPR (GPX4)",
        assay="MIX-Seq (pooled cell lines, SNP demultiplexing)",
        cell_context="Pools of 24 to 99 cancer cell lines", species="Homo sapiens",
        pub_perturbations="13 drugs plus one CRISPR target gene (GPX4)",
        pub_cells="Reported per experiment; no study-wide total stated",
        deposit="figshare 10.6084/m9.figshare.10298696", obtained_from="figshare",
        reference="McFarland et al. 2020, Nature Communications 11:4296",
        doi="10.1038/s41467-020-17440-w"),
    "Tahoe-100M": dict(
        modality="Small molecule",
        assay="Single-cell drug atlas (subset used here)",
        cell_context="50 cancer cell lines", species="Homo sapiens",
        pub_perturbations="1,100 small-molecule perturbations (379 distinct drugs)",
        pub_cells="~100 million transcriptomic profiles",
        deposit="Hugging Face tahoebio/Tahoe-100M and Arc Virtual Cell Atlas; "
                "no INSDC accession",
        obtained_from="Arc Virtual Cell Atlas",
        reference="Zhang et al. 2025, bioRxiv",
        doi="10.1101/2025.02.20.639398"),
}

# (dataset, field, confidence, url). One row per verified external claim.
P = "stated in primary source"
S = "stated in secondary source"
EXTERNAL_SOURCES = [
    ("Adamson", "cell context, modality, perturbations, accession", P,
     "https://pmc.ncbi.nlm.nih.gov/articles/PMC5315571/"),
    ("Adamson", "repository record GSE90546", P,
     "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90546"),
    ("Norman", "cell context, modality, 287 perturbations, ~110,000 cells, accession", P,
     "https://arcinstitute.org/work/science.aax4438.pdf"),
    ("Norman", "repository record GSE133344", P,
     "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133344"),
    ("Replogle", "cell lines, modality, 9,866 genes, >2.5M cells, SRA-only deposit", P,
     "https://europepmc.org/articles/PMC9380471"),
    ("Replogle", "BioProject PRJNA831566; no linked GEO series", P,
     "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA831566"),
    ("VCC", "H1 hESC, training split, hosting", P,
     "https://github.com/ArcInstitute/arc-virtual-cell-atlas/blob/main/virtual-cell-challenge/README.md"),
    ("VCC", "150 training targets and 183,097 cells (column sum of the official file)", P,
     "https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/"),
    ("VCC", "CRISPRi dual-guide modality", S,
     "https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge"),
    ("VCC", "venue, volume, pages, DOI", P,
     "https://api.crossref.org/works/10.1016/j.cell.2025.06.008"),
    ("sci-Plex", "three main-screen cell lines, 188 compounds, 649,340 cells", P,
     "https://pmc.ncbi.nlm.nih.gov/articles/PMC7289078/"),
    ("sci-Plex", "repository record GSE139944", P,
     "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139944"),
    ("McFarland", "24-99 cell-line pools, 13 drugs plus GPX4 CRISPR, no study-wide cell total", P,
     "https://pmc.ncbi.nlm.nih.gov/articles/PMC7453022/"),
    ("McFarland", "figshare deposit 10298696", P,
     "https://figshare.com/articles/dataset/MIX-seq_data/10298696"),
    ("Tahoe-100M", "50 cell lines, 1,100 perturbations, ~100M profiles", P,
     "https://www.biorxiv.org/content/10.1101/2025.02.20.639398v1.full"),
    ("Tahoe-100M", "hosting; no INSDC accession", P,
     "https://huggingface.co/datasets/tahoebio/Tahoe-100M"),
]

INTERNAL_SOURCES = [
    ("all allele datasets", "perturbations scored, cells scored, median depth",
     "data/pertresolve_bench.csv, via pertresolve.bench.PertResolveBench"),
    ("Adamson, Norman, Replogle, VCC, sci-Plex",
     "perturbations scored in the benchmark-ranking arm",
     "results/benchmark_resolution/summary.csv"),
    ("Adamson, Norman, Replogle and the four allele datasets",
     "perturbations scored in the rankability predictor",
     "results/rankability_predictor_honest.csv"),
    ("Adamson, Norman, Replogle", "perturbations in the prospective pilot arm",
     "results/fig6_derived/prospective_auc.csv"),
    ("Adamson, Norman, Replogle, VCC, McFarland, Tahoe-100M",
     "perturbations scored under the frozen criterion",
     "results/canonical/resolution_panel.csv"),
    ("Adamson, Norman, Replogle", "depth-matched atlas comparison membership",
     "results/fig6_derived/depthmatch_summary.csv"),
]

PANEL_KEY = {"adamson2016": "Adamson", "norman2019": "Norman", "replogle": "Replogle",
             "vcc_training": "VCC", "mcfarland": "McFarland", "tahoe100m": "Tahoe-100M"}

# Analysis arms. Membership is asserted here and CHECKED against the result
# files in build_roles(); a mismatch raises rather than being written out.
ARMS = {
    "Benchmark-ranking resolution (Fig. 5c)":
        ["TP53", "KRAS", "GATA1", "JAK1", "Adamson", "Norman", "Replogle", "VCC", "sciPlex"],
    "Rankability predictor, LODO (Fig. 6b,c)":
        ["TP53", "KRAS", "GATA1", "JAK1", "Adamson", "Norman", "Replogle"],
    "Pilot-to-evaluation validation (Fig. 6d)":
        ["Adamson", "Norman", "Replogle"],
    "Frozen packaged criterion (Fig. 6h)":
        ["Adamson", "Norman", "Replogle", "VCC", "McFarland", "Tahoe-100M"],
    "Native-vs-depth-matched atlas comparison (Fig. 6f)":
        ["Adamson", "Norman", "Replogle"],
}

DISPLAY = {"sciPlex": "sci-Plex"}
ORDER = ["TP53", "KRAS", "GATA1", "JAK1", "Adamson", "Norman", "Replogle",
         "VCC", "sciPlex", "McFarland", "Tahoe-100M"]


def _read(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"required result file missing: {path}")
    if path.suffix == ".json":
        return json.loads(path.read_text())
    return pd.read_csv(path)


def load_sources() -> dict:
    from pertresolve.bench import PertResolveBench
    bench_all = PertResolveBench.load(str(REPO / "data")).variants
    bench = bench_all[bench_all.variant.astype(str).str.upper() != "WT"].copy()
    return {
        "bench": bench,
        "bench_all": bench_all,
        "benchres": _read(RESULTS / "benchmark_resolution" / "summary.csv"),
        "panel": _read(RESULTS / "canonical" / "resolution_panel.csv"),
        "lodo": _read(RESULTS / "rankability_predictor_honest.csv"),
        "pilot": _read(RESULTS / "fig6_derived" / "prospective_auc.csv"),
        "depthmatch": _read(RESULTS / "fig6_derived" / "depthmatch_summary.csv"),
    }


def scored_counts(src: dict) -> dict[str, int]:
    """Largest perturbation count each external screen contributes to any arm."""
    benchres = src["benchres"].set_index("dataset")["n_perturbations"]
    lodo = src["lodo"][src["lodo"].feature_set == "effect_size"] \
        .set_index("held_out")["n_perturbations"]
    panel = src["panel"].copy()
    panel["name"] = panel.dataset_name.map(PANEL_KEY)
    panel = panel.set_index("name")["n_perturbations"]
    out = {}
    for name in ORDER:
        if name in ALLELE_STATIC:
            continue
        counts = [s[name] for s in (benchres, lodo, panel) if name in s.index]
        out[name] = int(max(counts))
    return out


def build_datasets(src: dict) -> pd.DataFrame:
    bench, bench_all = src["bench"], src["bench_all"]
    ext_n = scored_counts(src)
    rows = []
    for name in ORDER:
        if name in ALLELE_STATIC:
            st = ALLELE_STATIC[name]
            g = bench[bench.gene == name]
            g_all = bench_all[bench_all.gene == name]
            rows.append({
                "Dataset": name,
                "Layer": "Primary (allele-resolved)",
                "Perturbation modality": st["modality"],
                "Assay": st["assay"],
                "Cell context": st["cell_context"],
                "Species": st["species"],
                "Protein length (aa)": st["protein_len_aa"],
                "UniProt": st["uniprot"],
                "Perturbations as published": UNVERIFIED,
                "Cells as published": UNVERIFIED,
                "Perturbations scored here (n)": int(len(g)),
                "Cells scored here": int(g_all.n_cells.sum()),
                "Median cells per perturbation (here)": int(g.n_cells.median()),
                "Original data deposit": st["deposit"],
                "Obtained here from": st["obtained_from"],
                "Reference": st["reference"],
                "DOI": st["doi"],
            })
        else:
            st = EXTERNAL_STATIC[name]
            rows.append({
                "Dataset": DISPLAY.get(name, name),
                "Layer": "Broader perturbation screen",
                "Perturbation modality": st["modality"],
                "Assay": st["assay"],
                "Cell context": st["cell_context"],
                "Species": st["species"],
                "Protein length (aa)": "n/a",
                "UniProt": "n/a",
                "Perturbations as published": st["pub_perturbations"],
                "Cells as published": st["pub_cells"],
                "Perturbations scored here (n)": ext_n[name],
                "Cells scored here": UNVERIFIED,
                "Median cells per perturbation (here)": UNVERIFIED,
                "Original data deposit": st["deposit"],
                "Obtained here from": st["obtained_from"],
                "Reference": st["reference"],
                "DOI": st["doi"],
            })
    return pd.DataFrame(rows)


def build_roles(src: dict) -> pd.DataFrame:
    """Membership matrix, cross-checked against the files that define each arm."""
    observed = {
        "Benchmark-ranking resolution (Fig. 5c)":
            set(src["benchres"].dataset.str.replace("allele_", "", regex=False)),
        "Rankability predictor, LODO (Fig. 6b,c)":
            set(src["lodo"][src["lodo"].feature_set == "effect_size"].held_out),
        "Pilot-to-evaluation validation (Fig. 6d)": set(src["pilot"].dataset),
        "Frozen packaged criterion (Fig. 6h)":
            {PANEL_KEY[k] for k in src["panel"].dataset_name},
        "Native-vs-depth-matched atlas comparison (Fig. 6f)": set(src["depthmatch"].dataset),
    }
    for arm, declared in ARMS.items():
        if set(declared) != observed[arm]:
            raise AssertionError(
                f"arm membership disagrees with the result file for {arm!r}: "
                f"declared {sorted(declared)} vs observed {sorted(observed[arm])}")

    rows = []
    for ds in ORDER:
        row = {"Dataset": DISPLAY.get(ds, ds),
               "Layer": "Primary (allele-resolved)" if ds in ALLELE_STATIC
               else "Broader perturbation screen"}
        for arm in ARMS:
            row[arm] = "yes" if ds in ARMS[arm] else "-"
        row["Arms entered (n)"] = str(sum(ds in v for v in ARMS.values()))
        rows.append(row)
    total = {"Dataset": "Datasets in arm (n)", "Layer": "",
             **{arm: str(len(ARMS[arm])) for arm in ARMS}, "Arms entered (n)": ""}
    return pd.concat([pd.DataFrame(rows), pd.DataFrame([total])], ignore_index=True)


def build_definitions() -> pd.DataFrame:
    rows = [
        ("Layer", "Primary datasets are the allele-resolved screens in which the "
                  "measurement-resolution criterion is defined. Broader perturbation "
                  "screens test whether it holds beyond alleles."),
        ("Perturbations as published", "Scale reported by the original study. Verified "
                                       "from a primary source for the seven external "
                                       "screens only; see the Sources sheet."),
        ("Cells as published", "As above. Where a study reports only per-experiment "
                               "figures and no study-wide total, that is stated instead "
                               "of a total summed here."),
        ("Perturbations scored here (n)", "For allele datasets, the variant conditions in "
                                          "the benchmark. For external screens, the largest "
                                          "count across the arms it enters; the arms apply "
                                          "different inclusion thresholds and depths, so a "
                                          "single number would be misleading. Per-arm counts "
                                          "are in the manuscript Supplementary Tables."),
        ("Cells scored here", "Includes the dataset wild-type or control cells for the allele "
                              "datasets, matching the manuscript Methods."),
        ("Original data deposit", "Where the original study deposited the data, which is not "
                                  "always where this study obtained it."),
        ("Obtained here from", "The route actually used in this study."),
        (UNVERIFIED, "The value was not confirmed from a source that was retrieved. It is left "
                     "blank on purpose rather than filled from a secondary summary."),
        ("Results are not in this workbook",
         "Per-dataset measured quantities are in the manuscript Supplementary Information: "
         "Supplementary Table 5 (per-gene backbone and reproducibility ceiling), "
         "Supplementary Table 6 (signal and noise), Supplementary Table 7 (depth "
         "prescriptions) and Supplementary Table 9 (the pre-registered panel)."),
        ("Numbering", "Nature Portfolio numbers Supplementary Tables and Supplementary Data "
                      "as separate series. The Supplementary Information PDF holds "
                      "Supplementary Tables 1-10; this workbook is Supplementary Data 1."),
    ]
    return pd.DataFrame(rows, columns=["Term", "Definition"])


def build_sources() -> pd.DataFrame:
    rows = [("External metadata", d, f, c, u) for d, f, c, u in EXTERNAL_SOURCES]
    rows += [("Measured in this study", d, f, "computed from a committed result file", p)
             for d, f, p in INTERNAL_SOURCES]
    rows.append(("Primary allele datasets",
                 "assay, cell context, accession, protein length",
                 "stated in the manuscript Methods",
                 "manuscript/PertResolve_manuscript.pdf, Benchmark construction "
                 "and Data availability"))
    return pd.DataFrame(rows, columns=["Kind", "Dataset", "Fields", "Confidence", "Source"])


HEAD_FILL = PatternFill("solid", fgColor="1F3B4D")
HEAD_FONT = Font(bold=True, color="FFFFFF", size=10)
THIN = Side(style="thin", color="D0D7DC")
BAND = PatternFill("solid", fgColor="F4F6F7")


def write_workbook(sheets: dict[str, pd.DataFrame], out: Path) -> None:
    with pd.ExcelWriter(out, engine="openpyxl") as xl:
        for name, df in sheets.items():
            df.to_excel(xl, sheet_name=name, index=False)
            ws = xl.sheets[name]
            for cell in ws[1]:
                cell.fill, cell.font = HEAD_FILL, HEAD_FONT
                cell.alignment = Alignment(horizontal="center", vertical="center",
                                           wrap_text=True)
            ws.freeze_panes = "B2"
            ws.auto_filter.ref = ws.dimensions
            ws.row_dimensions[1].height = 32
            layer_col = [c.value for c in ws[1]].index("Layer") + 1 \
                if "Layer" in [c.value for c in ws[1]] else None
            for col in range(1, ws.max_column + 1):
                letter = get_column_letter(col)
                width = max(len(str(ws.cell(r, col).value or ""))
                            for r in range(1, ws.max_row + 1))
                ws.column_dimensions[letter].width = min(max(width + 2, 12), 44)
                for r in range(2, ws.max_row + 1):
                    c = ws.cell(r, col)
                    c.border = Border(bottom=THIN)
                    c.alignment = Alignment(vertical="top", wrap_text=True)
                    # Shade the primary allele rows so the two layers separate
                    # at a glance without adding a colour legend.
                    if layer_col and str(ws.cell(r, layer_col).value or "").startswith("Primary"):
                        c.fill = BAND
    print(f"wrote {out}")
    for name, df in sheets.items():
        print(f"  {name}: {df.shape[0]} rows x {df.shape[1]} cols")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, required=True, metavar="OUT_XLSX",
                    help="output .xlsx path; required, so no committed file is "
                         "overwritten by simply running the script")
    args = ap.parse_args()
    if args.out.suffix != ".xlsx":
        ap.error(f"--out must end in .xlsx, got {args.out}")
    reject_repo_results(args.out)

    src = load_sources()
    sheets = {
        "Datasets": build_datasets(src),
        "Analysis_roles": build_roles(src),
        "Column_definitions": build_definitions(),
        "Sources": build_sources(),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_workbook(sheets, args.out)


if __name__ == "__main__":
    main()
