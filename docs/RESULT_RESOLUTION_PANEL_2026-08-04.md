# Applying the frozen criterion to a wider panel: outcome against the pre-registration

Backs `results/canonical/resolution_panel.csv`. The criterion, the parameters, the datasets
and the predictions were fixed in `docs/PREREG_RESOLUTION_PANEL_v1.md` and committed at
`4c9a921`, before any dataset was measured. Nothing below was tuned after the fact and no
dataset was dropped.

## Result

Six perturbation screens, run through the installed package at `depth = 50`, `n_seeds = 8`,
`n_boot = 1000`, `seed = 0`, `n_components = 50`.

| dataset | perturbations | detectable | identifiable | replicate ceiling | P(order) | verdict |
|---|---|---|---|---|---|---|
| McFarland gene perturbation | 10 | 90.0% | 20.0% | 0.928 | 0.42 | detectable |
| Tahoe-100M demonstration | 7 | 71.4% | 42.9% | 0.929 | 0.11 | detectable |
| Norman 2019 | 195 | 68.7% | 0.5% | 0.898 | 0.997 | detectable |
| Replogle 2022 K562 essential | 414 | 51.2% | 1.4% | 0.865 | 0.96 | detectable |
| Adamson 2016 | 94 | 47.9% | 3.2% | 0.851 | 0.97 | not detectable |
| Virtual Cell Challenge training | 139 | 25.9% | 4.3% | 0.806 | 0.90 | not detectable |

**Not one of the six is benchmarkable** under the frozen criterion. The verdict ladder stops
at identification for every dataset in the panel.

## The four predictions

**1. Reproduction of the published ordering: failed in part.** The prediction was that the
detectable fraction would order Norman above Adamson above Replogle, matching the
manuscript's matched-50 rankable fractions of 88.8%, 61.5% and 45.1%. Norman is highest as
predicted, but Adamson and Replogle are reversed: 51.2% against 47.9%, a 3.3-point gap.

This is a real discrepancy and is reported as one. Three differences between the packaged
criterion and the manuscript's separate pipeline can produce it, and this run does not
distinguish them: the manuscript's pipeline chose Adamson's control automatically from a
candidate list while this run was given `Gal4-4(mod)_pBA582` explicitly; the inclusion
threshold here is `4 x depth = 200` cells against at least 100 there, which removed 1,643 of
Replogle's 2,058 perturbations and 21 of Adamson's; and the principal components are fit on
different cell subsets. The ordering of the extremes is preserved and the middle two are
within a few points of each other, so this does not overturn the framework, but the two
pipelines are not interchangeable and should not be quoted as if they were.

Adamson also sits almost exactly on the `DETECTABLE_FRACTION = 0.5` operating point at
47.9%, which is why its verdict reads "not detectable" while an earlier exploratory run at a
different control setting put it at 50.0% and one step higher. A verdict that turns on two
points either way is a reminder that these are declared operating points on continuous
quantities, not category boundaries.

**2. Detection without identification: confirmed on all six.** Every dataset has a far lower
identifiable than detectable fraction, and on the two largest screens the gap is extreme:
Norman is 68.7% detectable and 0.5% identifiable, Replogle 51.2% and 1.4%. The claim the
paper makes at allele resolution therefore extends to gene-level and drug-level screens.
Perturbations in these atlases can be told from the control and cannot be told from each
other at 50 cells per group.

**3. Set-level resolution survives per-perturbation failure: confirmed.** Replogle recovers
the full ordering of graded predictors with probability 0.96 while 1.4% of its perturbations
are identifiable; Norman 0.997 at 0.5%; the Virtual Cell Challenge set 0.90 at 4.3%. A large
screen aggregates many individually unresolved perturbations into a benchmark that still
orders predictors, which is the dissociation the manuscript reports, reproduced here by an
independent implementation on datasets it was not built on.

The two small panels behave in the opposite direction and confirm the same mechanism from
the other side: Tahoe with 7 perturbations recovers the ordering with probability 0.11 and
McFarland with 10 with probability 0.42, despite having the highest identifiable fractions
and the highest ceilings in the panel. Individually resolvable perturbations do not make a
benchmark if there are too few of them.

**4. No benchmarkable verdict with near-zero identification: confirmed.** There were no such
verdicts, because there were no benchmarkable verdicts at all.

## What this does and does not establish

It establishes that an implementation frozen before the run, applied without per-dataset
tuning to six screens, three of which took no part in developing it, reproduces the
detection-without-identification ordering and the per-perturbation versus set-level
dissociation, and that no screen in the panel supports identification of individual
perturbations at this depth.

It does not establish agreement with the manuscript's own pipeline at the level of
individual fractions, as prediction 1 shows. It also rests on six datasets, below the ten to
fifteen the plan anticipated, because six is the number of usable perturbation screens on
this machine; two of the six carry fewer than a dozen perturbations and their fractions are
correspondingly uncertain.
