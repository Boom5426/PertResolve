# Pre-registration: applying the frozen resolution criterion to a wider panel

Written before the panel is run. Its purpose is to make the run capable of contradicting
the framework rather than illustrating it, by fixing what will be measured, on what, and
what result would count against the claim, while none of those answers is known.

The four allele datasets and the five external screens already in the manuscript were used
to develop this criterion. They cannot test it. This applies the frozen criterion to
datasets that took no part in building it.

## What is frozen

The criterion is `alleleperturb.resolution.resolution_report` at commit `dea8225`, run
through the installed package with no per-dataset tuning. Specifically:

- four disjoint groups of `depth` cells per perturbation, cut from one shuffle;
- detection: split-half energy distance on cells, a perturbation counts as detectable when
  its distance from the control exceeds the full 95% spread of its own replicate distance;
- identification: the cross-fitted nearest-competitor separation must exceed its own 95%
  spread across cell splits, the same margin detection carries;
- the replicate ceiling: tie-aware cosine discrimination of a second measurement of the
  same perturbation;
- ordering recovery: graded predictors interpolating toward another perturbation's profile,
  bootstrapped over perturbations;
- verdict thresholds `DETECTABLE_FRACTION = 0.5`, `IDENTIFIABLE_FRACTION = 0.5`,
  `BENCHMARKABLE_CEILING = 0.65`.

Run parameters, fixed now: `depth = 50`, `n_seeds = 8`, `n_boot = 1000`, `seed = 0`,
`n_components = 50`. No parameter will be changed after seeing any result. If a dataset
fails to run, the failure is reported; it is not rescued by changing a setting.

## The panel

Every perturbation dataset on the compute server that carries a perturbation column, a
usable control group and at least five perturbations with 200 cells. Surveyed before any
measurement:

| dataset | cells | perturbation column | control | perturbations with 200 cells |
|---|---|---|---|---|
| Norman 2019 | 111,445 | `perturbation` | `control` | 196 |
| Replogle 2022 K562 essential | 310,385 | `perturbation` | `control` | 415 |
| Adamson 2016 | 65,337 | `perturbation` | `Gal4-4(mod)_pBA582` | 95 |
| Virtual Cell Challenge training | 221,273 | `target_gene` | `non-targeting` | 140 |
| McFarland gene perturbation | 5,500 | `perturbation` | `control` | 11 |
| Tahoe-100M demonstration subset | 80,000 | `perturbation` | `DMSO` | 8 |

Excluded, with the reason, before running: `cigs_hek293t_with_controls` (10,726
perturbations, only 2 with 200 cells); `Parse_10M_PBMC_cytokines` (two treatment groups, not
a screen); `perturb_processed` (64 cells); the sci-Plex data on this machine is unprocessed
raw text and was not converted for this run.

**Three of these six are in the manuscript already** (Norman, Replogle, Adamson) and are
included deliberately, as a reproduction check: the packaged criterion should place them
where the manuscript's separate pipeline did. The Virtual Cell Challenge set is in the
manuscript's benchmark-resolution figure but has not been through this criterion. McFarland
and Tahoe are new to this work entirely.

Six is below the ten to fifteen the plan anticipated. That is the number of usable
perturbation screens present on this machine, and the panel is not padded with datasets that
cannot support the measurement. Extending it requires downloads that are outside this run.

## Predictions, recorded now

1. **Reproduction.** Norman, Replogle and Adamson will rank in the same order on the
   detectable fraction as in the manuscript's matched-50 analysis, which put un-rankable
   fractions at 11.2%, 54.9% and 38.5% respectively, so detectable should be highest for
   Norman and lowest for Replogle. A different ordering would mean the packaged criterion is
   not measuring what the manuscript's pipeline measured.

2. **Detection without identification.** Most datasets will clear detection for a
   substantial share of perturbations and clear identification for far fewer. This is the
   claim the paper makes at allele resolution, extended to gene and drug level. A dataset
   whose identifiable fraction matches or exceeds its detectable fraction would contradict
   the ordering the whole framework rests on.

3. **Set-level resolution survives per-perturbation failure.** Ordering recovery will be
   high on the large screens even where the identifiable fraction is low, reproducing the
   dissociation between per-perturbation and set-level resolution.

4. **No dataset will be reported benchmarkable at allele-style resolution while having a
   near-zero identifiable fraction**, because the verdict ladder checks identification
   before the ceiling.

## What would count against the framework

- Prediction 1 failing: the criterion does not reproduce the published ordering.
- Prediction 2 failing on more than one dataset: identification is not the harder question,
  and the ladder is wrong.
- The verdict disagreeing with the underlying numbers on any dataset, for instance a
  benchmarkable verdict with an identifiable fraction near zero.
- Any dataset on which the criterion cannot be computed for a reason the code does not
  report clearly.

## What will be reported either way

Every dataset in the table, including any that fail, with the verdict, all three fractions,
the ceiling, both signal statistics, and the number of perturbations excluded for want of
cells. No dataset will be dropped after the fact, and no threshold will be moved.
