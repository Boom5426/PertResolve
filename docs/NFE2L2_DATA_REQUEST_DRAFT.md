# NFE2L2 Saturation-seq: data request draft

**Not sent.** Draft only, for the author to review, edit and send.

## Why a request is needed

Searched and found nothing public:

| source | query | result |
|---|---|---|
| GEO (eutils, `db=gds`) | "Saturation-seq", "NFE2L2 variant single cell", "saturation genome editing single-cell RNA" | no matching series |
| ENA portal API | `study_title="*Saturation-seq*"`, `"*NFE2L2*"` | no matching study |
| BioStudies | "Saturation-seq", "Saturation-seq NFE2L2" | 0 hits |
| Europe PMC data links | `PPR/PPR1271599/datalinks` | hitCount 0 |
| bioRxiv full text | `10.64898/2026.06.30.735631` | Cloudflare 403, statement unreadable from here |
| **GitHub `StatOVarI-lab/Saturation-seq`** | 65 files, README, all shell and R scripts | **code only, no accession, no data directory, no download URL** |

The repository is public, last pushed 2026-06-26, described as "Code for analysis
of single-cell-based saturation genome editing with transcriptomic readout". Its
README lists `run_DNA_processing.sh` and `RNA_analysis.sh` as the entry points and
says nothing about where the data lives. Grepping every script for `GSE`, `PRJ`,
`EGA`, `E-MTAB`, `zenodo`, `figshare`, `http`, `ftp` returns nothing.

So the processed data is almost certainly not yet public, and the bioRxiv Data
availability statement should still be read directly before sending, in case it
names an accession under embargo or an "available on request" clause to cite.

## Recipients

- **David J. Adams**, corresponding author, Wellcome Sanger Institute (per the
  bioRxiv record).
- **Magdalena E. Strauss**, first author; `StatOVarI-lab` appears to be their
  group, so they are the likely practical contact for the processed matrices.
- Consider copying **Andrew J. Waters** and **Sarah Cooper** (Sanger), who are on
  the author list and on the experimental side.

Addresses are not recorded here; take them from the preprint's corresponding
author footnote.

## Draft

> **Subject:** Data request: Saturation-seq NFE2L2 processed matrices for a
> variant-effect prediction benchmark
>
> Dear Dr Adams and Dr Strauss,
>
> I read your preprint "Saturation-seq integrates single-cell saturation genome
> editing and RNA-seq to quantify NFE2L2 (NRF2) variant effects"
> (doi:10.64898/2026.06.30.735631) with great interest. Installing 230 variants at
> the endogenous locus in a barcoded haploid background and linking each edit to a
> single-cell transcriptome is exactly the design that the field has been missing.
>
> I am completing a benchmark study of whether current single-cell perturbation
> models can predict the transcriptional consequences of individual protein-coding
> variants. Our finding so far is a negative one: across four allele-resolved
> datasets, models sit at chance, and in three of the four the measurement itself
> cannot reproducibly distinguish sibling alleles, so the benchmark cannot reward a
> correct answer even in principle. Only one dataset clears that bar, which leaves
> the conclusion resting on a single positive example. A second gene with many
> alleles in a uniform background would let us test whether the model failure
> generalises, and your NFE2L2 panel is the closest match we have found.
>
> I could not locate an accession for the data in GEO, ENA, BioStudies or Europe
> PMC, and the analysis repository at github.com/StatOVarI-lab/Saturation-seq
> contains code without a data pointer, so I wanted to ask directly rather than
> assume. Would you be willing to share, at whatever stage suits you:
>
> 1. the processed single-cell RNA count matrix (raw or normalised counts, cells by
>    genes);
> 2. the per-cell variant assignment from the amplicon arm, including the
>    unassigned and multiply-assigned cells rather than only the confident calls;
> 3. the variant annotation table: HGVS or protein-level identifier, edit class,
>    and any quality or confidence flag you use for filtering;
> 4. replicate and batch metadata: which cells belong to which biological
>    replicate, transduction pool, plate or lane.
>
> Item 4 matters more than it may appear. Our analysis turns on separating
> reproducible allele-specific signal from batch structure, and without replicate
> and batch labels an apparent allele effect cannot be distinguished from a
> position or pool effect.
>
> On terms: we are happy to work under whatever arrangement you prefer, including
> an embargo until your paper is published, a data transfer or use agreement, or
> restricting our use to the specific analyses described above. We would of course
> cite the preprint, and we would be glad to share our analysis code and results
> with you before any submission, and to discuss authorship or acknowledgement if
> you would find that appropriate. If the data is already scheduled for public
> release, simply pointing us at the accession and timeline would be just as
> useful.
>
> Thank you for considering this, and congratulations on a very nice piece of work.
>
> With best regards,
>
> Bo Li
> [affiliation, position]
> [email]

## Before sending, please check

- **Read the preprint's own Data availability statement.** If it names an
  accession or says "available upon request", quote it back in the email; if it
  promises release on publication, the request becomes a timeline question rather
  than an access question.
- Decide in advance what you are willing to offer: co-authorship, acknowledgement,
  or citation only. The draft leaves this open, which is fine for a first contact
  but should not stay open in a follow-up.
- Confirm the affiliation line and add your institutional email; a request from an
  institutional address is answered far more often than one from a personal one.
- If there is no reply in about two weeks, a short follow-up to the first author
  alone is usually more effective than another message to the whole list.

## What a positive reply would unlock

230 endogenous-locus variants in one gene and one background, with paired
genotype and transcriptome per cell, is a substantially larger and cleaner allelic
series than anything currently in the benchmark: the existing JAK1 arm has 26
variants and the TP63 bulk arm has 18. It is the only realistic route to testing
whether the model-limited finding generalises beyond a single gene.
