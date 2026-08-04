# The residual scoring axis, and why 0.5 is the wrong reference on it

Backs `results/canonical/residual_axis_{models,ceiling,per_gene}.csv`, produced by
`scripts/analysis/residual_axis.py`. Adds a second scoring axis; changes no existing
definition and retrains no model.

## Why a second axis

Most of a variant's measured response is the programme its gene shares with every sibling,
so a prediction that reproduces only that programme is right about the gene and silent about
the allele. The Discussion already says predictors "should be judged on whether they recover
allele-specific residual structure, not the gene-shared response that dominates
direction-based scores", and nothing in the paper scored that way. This does, for all 25
methods, on the same held-out variants and candidate sets, from the same stored predictions.

## The gene mean now comes from training variants only

`results/reviewer_controls/residual_decomp_*.csv` centred on a mean that included the
variants being scored: over every variant on the measurement side, and over the test
variants alone on the model side. Recomputed with a training-only mean, the replicate
ceiling on residuals for JAK1 falls from 0.892 to **0.839**, so the earlier figure was
optimistic, as self-inclusion in the mean predicts.

## The finding that matters: the residual axis does not have chance at 0.5

`WT-null` predicts a zero vector. On the residual axis that becomes minus the gene mean, one
constant shared by every held-out variant, carrying no information whatever. It scores
**0.524**, and its permutation null is **0.524** exactly, because permuting a constant
changes nothing.

The reason is structural. The mean removed is built from training variants, so their
residuals are shrunk toward the origin by the 1/n share of themselves they contain, while
held-out variants' residuals are not. A constant query therefore ranks held-out variants
systematically differently from training ones, and the axis has an offset that no model
produced.

Anything read against 0.5 on this axis is therefore misread. Taken at face value, seven of
twenty-five methods have residual-PDS intervals excluding 0.5. Against their own permutation
nulls, the picture is different.

## Against their own permutation nulls

Predictions are permuted across the held-out variants of a gene, preserving every marginal
property of predictions and truths and destroying only which goes with which.

| method | residual-PDS | its null | excess | P |
|---|---|---|---|---|
| MLP-theta | 0.538 | 0.513 | 0.025 | 0.020 |
| Ridge-theta | 0.525 | 0.502 | 0.023 | 0.020 |
| Biolord | 0.521 | 0.504 | 0.017 | 0.020 |
| GBoost-theta | 0.528 | 0.512 | 0.016 | 0.039 |
| RF-theta | 0.527 | 0.512 | 0.015 | 0.039 |
| scGen | 0.528 | 0.517 | 0.011 | 0.020 |

Six of twenty-five reach P < 0.05. **None is established**, for three reasons that should be
stated together rather than one at a time. The largest excess is 0.025, two and a half
points of a score whose range is half the unit interval. Fifty permutations put a floor of
1/51 = 0.0196 on the attainable P value, and four of the six sit exactly on that floor, so
they are "best of fifty" rather than resolved. And twenty-five methods were tested: a
Bonferroni threshold would be 0.002, which fifty permutations cannot produce at all.

The honest reading is that the residual axis leaves open the possibility of a small
allele-specific signal in the theta-feature heads and two external models, and settles
nothing. Every one of the six is a theta-based head or an external model and none is an
ESM-based head, which is suggestive and no more.

## The null result is unchanged, and the model gap widens

| gene | ceiling, PDS | ceiling, residual-PDS |
|---|---|---|
| TP53 | 0.473 (0.445 to 0.499) | 0.509 (0.488 to 0.531) |
| KRAS | 0.492 (0.450 to 0.537) | 0.518 (0.487 to 0.551) |
| GATA1 | 0.581 (0.560 to 0.602) | 0.576 (0.557 to 0.596) |
| JAK1 | 0.779 (0.730 to 0.826) | **0.839 (0.792 to 0.883)** |

Read against the 0.524 offset the axis carries, the TP53 and KRAS residual ceilings sit at
or below it: removing the shared programme exposes nothing there, which is what a
measurement-limited dataset should do. GATA1 is modestly above. JAK1 rises from 0.779 to
0.839 while the best model gains at most 0.025, so the residual axis **widens** the distance
between what the measurement permits and what the models reach. That is the direction that
matters: an axis that made the models look better would be the one to distrust.

## What must accompany any use of this axis

State the permutation null next to every number. A residual-PDS of 0.53 is above 0.5 and at
its null; the two readings differ in what they mean and only one of them is about a model.
