# Candidate second model-limited system: verification audit (2026-08-03)

Scope: verify four externally proposed candidate datasets and, where the data is
public, measure whether the measurement resolves allele identity. Nothing here is
a manuscript claim yet. Every number below was produced by a script in
`scripts/dataset_screen/` from data downloaded from the accession named beside it.

Standard applied, from the existing framework: a system is *model-limited* only if
**replicate/split-half identification is above chance while models sit at chance**.
This audit settles the first half. The second half (do ESM / biophysical features
predict the allele residual?) is a separate run and is NOT answered here.

---

## 1. Verification status of the four proposals

| # | Candidate | Exists? | Public data? | Verdict |
|---|---|---|---|---|
| 1 | TP63 SCRAMseq (GSE311874 family) | yes | partly | **usable now, but not the part that was proposed** |
| 2 | NFE2L2 Saturation-seq (Strauss et al.) | yes | **not found** | blocked on data availability |
| 3 | PAX5 CRAFT-seq (Baglaenko et al.) | yes | yes (mixed) | **wrong variant class, reject** |
| 4 | SDR-seq (GSE268646) | yes | yes | **targeted RNA panel, reject** |

Two errors in the source material worth recording:

- GSE311882 is **`scRNA_MYOD1`**, not `scRNA_P63noconv`. The P63 non-converted
  series is GSE311881.
- The study is **eight** subseries (GSE311874, 311876-311882), not three. The
  source listed only the umbrella plus two, and missed **GSE311877**, which is
  the single most useful piece of the whole deposit.

The `10.64898` bioRxiv DOI prefix, which looked like a fabrication signal, is
genuine: bioRxiv's current prefix. Confirmed against Crossref and the bioRxiv API.

---

## 2. TP63 SCRAMseq: what is actually in the deposit

Eight subseries, all released 2026-07-10, submitted 2025-12-01.

| Accession | Content | n | Processed files |
|---|---|---|---|
| GSE311874 | MITE_DBD amplicon (MAVE) | 8 | RAW.tar, 120 KB |
| GSE311876 | MITE_SAM amplicon (MAVE) | 8 | RAW.tar, 100 KB |
| GSE311877 | **RNAseq_mutants, arrayed bulk** | **80** | **RAW.tar, 9.0 MB** |
| GSE311878 | RNAseq_MYOD1 bulk | 8 | - |
| GSE311879 | RNAseq_P63 bulk | 8 | - |
| GSE311880 | scRNA_P63conv, Smart-seq3Express | 2 | RAW.tar, 10.9 GB |
| GSE311881 | scRNA_P63noconv | 4 | RAW.tar, 4.8 GB |
| GSE311882 | scRNA_MYOD1 | 4 | RAW.tar, 4.9 GB |

### 2.1 The scRNA arm cannot be used as released

Exact cell counts, read from the released barcode files:

| Sample | Pool | Cells |
|---|---|---|
| GSM9432943 | TP63 DBD | 18,432 |
| GSM9432944 | TP63 SAM | 18,416 |

18,432 = 48 x 384, i.e. plate-based Smart-seq3Express.

**Blocker.** The processed supplementary files are exactly three per sample:
`*.dgecounts.rds` (zUMIs), `*.gene_names.txt.gz`, `*kept_barcodes.txt.gz`.
There is **no per-cell variant assignment table**. Variant identity has to be
called from full-length reads over the TP63 transgene, which means reprocessing
the raw FASTQs from SRA (BioProject **PRJNA1371870**) against a custom construct
reference. That is the entire cost of using this arm, and it is not small.

**How many variants would survive if that reprocessing worked.** The MAVE
amplicon tables give the pool composition, so the cells-per-variant distribution
can be estimated rather than guessed:

| Pool | Variants | Residues | Range | Cells | >=5 cells | >=10 | >=20 | median cells |
|---|---|---|---|---|---|---|---|---|
| DBD (converted) | 1,179 | 60 | 217-276 | 18,432 | 744 | 475 | 136 | 7.8 |
| DBD (non-conv.) | 1,179 | 60 | 217-276 | 18,432 | 904 | 562 | 139 | 9.6 |
| SAM (converted) | 1,074 | 56 | 470-525 | 18,416 | 755 | 537 | 250 | 10.1 |

Near-saturation: ~19-20 substitutions per residue. Total 2,253 variants, matching
the "~2,300 missense variants" in the GEO summary. So the optimistic reading of
the source material was directionally right about scale, but the realistic
per-variant depth is **~8-10 cells median, ~500 variants at >=10 cells**, not
thousands of well-powered alleles. Top 10% of variants hold 46-56% of the pool.

Open question flagged, not resolved: the MAVE residue ranges (217-276, 470-525)
do not line up with the arrayed mutant positions (279-311, 514-556) under a
single offset. Likely a TAp63/dNp63 isoform numbering difference. Must be settled
against the paper before any cross-referencing of the two arms.

### 2.2 GSE311877 is the usable piece, and it was missed

80 bulk RNA-seq libraries = **18 TP63 missense alleles + WT + GFP, 4 biological
replicates each**, one gene, one uniform fibroblast-to-iKC conversion assay,
33,694 Ensembl genes, library size 1.6-5.0 M (median 2.4 M). Processed counts,
9 MB, no reprocessing needed.

Alleles, two regions:
- DBD-side: R279Q, R279S, G293H, G297L, R304Q, R304T, G310Q, R311I
- SAM/TID-side: L514D, L514F, C522D, C522G, T527I, L531E, L531R, L548N, I550D, H556A

Five **same-residue sibling pairs**: R279Q/S, R304Q/T, L514D/F, C522D/G, L531E/R.
This is the hardest identification test the design supports, and it exists here.

---

## 3. Does the measurement resolve allele identity? (GSE311877)

`scripts/dataset_screen/audit_gse311877.py`. Split-half over allele deltas from
WT, averaged over all three balanced 2v2 replicate splits. The WT reference is
split alongside the mutants (half A uses WT reps of half A), so no shared WT term
is injected into both halves. Gene selection is by mean expression, which uses no
allele labels. 12,734 genes kept of 33,694.

```
chance PDS = 0.500          chance top-1 = 0.056  (18 alleles)

split (1,2) v (3,4)   PDS 0.948   top-1 0.611
split (1,3) v (2,4)   PDS 0.967   top-1 0.778
split (1,4) v (2,3)   PDS 0.977   top-1 0.833
mean                  PDS 0.964   top-1 0.741
averaged matrix       PDS 0.977   permutation null 0.501, p95 0.611, p = 0.0005
```

Single-replicate delta correlations: within-allele r = 0.775 (n=108),
between-allele r = 0.597 (n=2,448). The gap is real but modest, which is the
warning sign that most of the delta is shared across alleles.

### 3.1 Is that allele resolution, or just severity ranking?

`scripts/dataset_screen/residual_gse311877.py`. This is the test that matters,
because a single shared failure axis with allele-specific magnitude would produce
a high PDS while carrying no allele-specific direction at all.

```
identify an allele from ||delta|| alone (ONE scalar):
                              PDS 0.836   top-1 0.259

variance of the allele-delta space:
    PC1 0.467  PC2 0.172  PC3 0.077  PC4 0.063  PC5 0.036  PC6 0.028

split-half identification after projecting out the top-k shared PCs
(PCs estimated on half A only, then applied to both halves):
    k        PDS     top-1    perm null    p
    0       0.964    0.741      0.501    0.0005
    1       0.806    0.519      0.502    0.0005
    2       0.817    0.481      0.497    0.0005
    3       0.813    0.481      0.499    0.0005
    5       0.727    0.278      0.501    0.0005
```

Reading: severity explains most of the headline number (a single scalar already
scores 0.836 of the 0.964). But **allele-specific residual structure survives
removal of the top five shared components**: PDS 0.727, top-1 0.278 against a
0.056 chance rate, permutation p = 0.0005.

That is materially different from the current AllelePerturb systems, where the
residual decomposition put residual Pearson near zero and the TP53/KRAS
residual-oracle at chance. Same-residue siblings all self-match above their
sibling, but by thin margins (R304Q self 0.735 vs cross 0.713).

### 3.2 The design confound, and why it does not appear to bite

`scripts/dataset_screen/confound_gse311877.py`. Layout recovered from the GEO
filenames: 12 alleles sit on one 96-well plate, **one allele per (column, row-block)**,
replicates running down the rows. The remaining 6 alleles plus WT and GFP have no
well id and form a separate batch. So **allele identity is aliased to well position
by construction**, and no allele is replicated across positions. That cannot be
fixed from the released metadata.

It can, however, be tested. Six columns each host two different alleles. If a
column artifact were carrying the signal, those pairs should be unusually similar:

```
                          same-column pairs   all other pairs   diff      p
raw delta space               r 0.8408          r 0.8216      +0.0191   0.57
residual space (top-5 PC out) r 0.8363          r 0.8396      -0.0033   0.83

row-block structure, residual space:
    rows A-D (n=7)   within 0.8424   across 0.8291
    rows E-H (n=5)   within 0.8426   across 0.8299
    no well id (n=6) within 0.7561   across 0.8225
```

No detectable column effect in either space. Block effects are +0.013 on a 0.83
baseline, and the no-well batch is *less* internally similar than across, i.e.
anti-clustered. **The confound is real in the design but leaves no measurable
footprint in the data.** State it as a limitation; do not claim it is absent.

Second design note: WT and GFP are in the no-well batch, so every delta for the
12 plated alleles carries a constant cross-batch offset. It is identical across
those alleles and therefore cancels under the column-centering used here, but it
inflates absolute effect sizes versus WT.

---

## 4. The other three candidates

**NFE2L2 Saturation-seq** - Strauss, Waters, Robertson, ..., Adams (Wellcome
Sanger), bioRxiv `10.64898/2026.06.30.735631`, posted 2026-07-04. Real; confirmed
via the bioRxiv API and Crossref. 230 variants in the N-terminal region of
NFE2L2, installed at the endogenous locus in a barcoded haploid line, read out by
single-cell amplicon plus transcriptome. Conceptually the closest match to what
this project needs.

**No public accession could be located.** Searched GEO (eutils), ENA (portal
API), BioStudies, and the Europe PMC data-links service: all negative. bioRxiv
itself is Cloudflare-blocked from this environment, so the Data availability
statement could not be read. Next step is manual: open the PDF, read the
statement, and if there is no accession, contact the corresponding author.

**PAX5 CRAFT-seq** - Baglaenko et al., *Nature* 646:117-125 (2025),
PMC12488502. Data availability, verbatim: count matrices at Zenodo
`10.5281/zenodo.15935857`, raw sequencing at dbGaP `phs004042.v1.p1` (controlled
access). Code at `github.com/immunogenomics/craft-seq`.

Reject. The paper's subject is gene disruption, regulatory-region deletions and
**non-coding SNPs** (the headline result is IL2RA `rs61839660` in primary T
cells); the PAX5 component is multiplexed editing of two regions, not a missense
allelic series. AllelePerturb is a protein-coding missense benchmark, so this is
the wrong variant class, independent of the dbGaP access barrier.

**SDR-seq** - GSE268646, "Functional phenotyping of genomic variants using joint
multiomic single-cell DNA-RNA sequencing". 17 samples, human and mouse, released
2025-06-17, 463 MB of MTX/TSV. Reject: the RNA side is a **targeted panel** (up
to 480 loci and their associated genes), not the transcriptome-wide response
AllelePerturb evaluates, and the modules are CRISPRi/PE/BE method demonstrations
rather than a same-gene allelic series.

---

## 5. Where this leaves the plan

What is settled: GSE311877 gives a same-gene, uniform-background, 4-replicate,
18-allele coding series where **the measurement resolves allele identity well
above chance, and does so beyond the shared severity axis**. That is a genuine
second measurement-resolved system, obtainable today, at 9 MB.

What is not settled, in order of importance:

1. **Are models at chance on it?** Untested. This is the half that decides
   whether it is *model-limited* rather than merely resolved. Needs the existing
   harness plus ESM / biophysical features over the 18 alleles.
2. **18 alleles is thin for held-out prediction.** Leave-one-variant-out over 18
   points is underpowered for a feature-to-response regression; the honest use of
   GSE311877 may be as a *measurement-ceiling* system rather than a prediction
   benchmark. Decide this before building a claim on it.
3. **The scRNA arm needs SRA reprocessing** to yield ~500 variants at >=10 cells.
   Large effort, uncertain yield, and the only route to a large allelic series
   from this study.
4. **NFE2L2 data availability** is a five-minute manual check that could change
   the ranking entirely.

Recommended ordering: 1 before 3. If models are not at chance on the 18 alleles,
the expensive scRNA reprocessing is not worth starting.

---

## Reproduction

From the repository root:

```bash
mkdir -p data_external/GSE311877 && cd data_external/GSE311877
curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311877/suppl/GSE311877_RAW.tar
tar -xf GSE311877_RAW.tar && rm GSE311877_RAW.tar && cd ../..

python3 scripts/dataset_screen/audit_gse311877.py       # section 3
python3 scripts/dataset_screen/residual_gse311877.py    # section 3.1
python3 scripts/dataset_screen/confound_gse311877.py    # section 3.2
```

`data_external/` is gitignored; the downloaded counts are never committed.
Override the input and output locations with `GSE311877_DIR` and
`DATASET_SCREEN_OUT` if needed. Seed 0 throughout; 2,000 permutations for the
identification nulls, 20,000 for the column-confound null. Outputs land in
`results/dataset_screen/`. Requires numpy and pandas only. All three scripts were
run to completion from a clean checkout of these paths; the numbers in sections 3,
3.1 and 3.2 are their actual stdout.
