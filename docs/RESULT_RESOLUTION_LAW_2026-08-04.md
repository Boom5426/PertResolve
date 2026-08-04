# The measurement quantity that governs discrimination, and what it prescribes

Backs `results/canonical/resolution_sweep.csv`, `resolution_law.json` and
`resolution_law_prescriptions.csv`, produced by `scripts/analysis/resolution_sweep.py`
and `scripts/analysis/resolution_law.py`. Continues
`docs/RESULT_COLLAPSE_REFUTED_2026-08-04.md`, which established that the quantity Methods
names does not govern discrimination. This identifies the one that does and inverts it
into a depth requirement.

## Design

1,452 points. Twelve are the four allele datasets at three depths, measured with no
construction. The other 1,440 add a synthetic separation of controlled size **and
controlled geometry** to TP53 and KRAS cells, whose own between-variant signal is measured
at `|rho2| < 0.005` with every interval covering zero. Rescaling measured profiles cannot
serve as a dial, since it scales signal and noise together and leaves their ratio fixed;
the dial has to act on the signal alone.

Each variant's cells are cut into five disjoint groups from one shuffle: a direction block
that nothing else touches, a query and a truth for the ceiling, and two that estimate the
signal. The two axes of every point therefore come from different cells. An earlier draft
drew the direction block from its own shuffle, leaving it overlapping the scoring cells; a
direction correlated with the noise it is added to would have produced the relationship by
construction rather than measuring it.

## Which summary of the measurement governs discrimination

| axis | usable points | Spearman with the ceiling |
|---|---|---|
| `rho2`, the mean squared separation over pairs, what Methods names | 1,452 | 0.767 |
| **`rho2_nn_median`, the median nearest-competitor separation** | 1,452 | **0.961** |
| `rho2_nn_geomean` | 329 | 0.818 |
| `frac_above_noise` | 1,452 | 0.950 |

On the twelve real points alone `rho2` scores 0.965 and the nearest-competitor median
0.881, so **four datasets cannot distinguish the two**. That is why the assertion survived
until it was tested against configurations built to separate them.

The geometric mean is unusable and is reported as undefined rather than floored: 1,123 of
1,452 points contain at least one perturbation whose separation debiases below zero, and a
floored geometric mean is set by the floor and by how many values hit it while still
printing as a number.

The above-noise fraction is close behind the chosen axis but is computed per cell split and
then averaged, never from separations averaged first. Thresholding an average shrinks the
noise around a residual bias rather than the bias itself, which on data with exactly zero
separation reads 0.59 instead of 0.52; the table was regenerated after that was fixed, and
the fraction's agreement with the ceiling moved from 0.923 to 0.950 while the chosen axis
and every conclusion below were unaffected.

## What the calibration can and cannot be inverted

Only the ceiling is monotone in the chosen axis (Spearman +0.838 on the lower half of the
axis, +0.918 on the upper). The two ordering-recovery outcomes **turn over**:

| outcome | lower half | upper half | peak |
|---|---|---|---|
| `ceiling_pds` | +0.838 | +0.918 | monotone |
| `p_correct` | +0.780 | **-0.280** | `rho2_nn_median` about 0.043, ceiling about 0.705 |
| `p_winner` | +0.757 | **-0.248** | same |

The turnover is real rather than an artefact. At the floor no predictor is distinguishable
from another; at high resolution the best two graded predictors both score near the
ceiling and are again indistinguishable, so the probability of recovering their full
ordering falls. A calculator fitted to a monotone form would therefore be confidently wrong
at the top, which is why monotonicity is tested before anything is inverted.

This bears on the manuscript. Results currently state that above a ceiling of 0.65 both
recovery measures "approached 1.0" (Fig. 5f, 5g). That describes the rising limb, and the
peak sits at a ceiling near 0.705, just above it. The claim is not contradicted within the
range the four datasets occupy, but it should not be extended upward.

## The design calculator, and the three cases where it refuses

The nearest-competitor separation belongs to the perturbations and does not change with
depth; the noise does, as `eta2 ~ 1/m`, confirmed empirically in Phase 1. So resolution
scales linearly in depth from a pilot, and the depth reaching a target ceiling is read off
the calibration by inversion.

It returns nothing rather than a number when the calibration never reaches the target at
that candidate-set size, when the pilot estimate is not positive, or when the pilot's lower
confidence bound does not clear zero. The third gate is the one that matters:

| dataset | pilot depth | `rho2_nn_median` (95% lower bound) | cells per group for a ceiling of 0.7 | for 0.9 |
|---|---|---|---|---|
| TP53 | 25, 50, 100 | -0.0026 (-0.0058), 0.0005 (-0.0079), 0.0019 (-0.0033) | refused | refused |
| KRAS | 25, 50, 100 | 0.0004 (-0.0056), 0.0010 (-0.0063), -0.0073 (-0.0133) | refused | refused |
| GATA1 | 25 | 0.0157 (0.0027) | 62 | 232 |
| GATA1 | 50, 100 | 0.0010 (-0.0108), 0.0051 (-0.0141) | refused | refused |
| JAK1 | 25, 50, 100 | 0.075 (0.049), 0.127 (0.091), 0.167 (0.142) | **10, 12, 19** | **55, 65, 99** |

Without that gate the calculator issued 2,000 to 15,000 cells per perturbation for TP53 and
KRAS, extrapolated from estimates of 0.0005 whose intervals cover zero. That reads as a
plan and is an extrapolation from noise. No finite depth rescues a separation that is not
there, and the honest output is a refusal.

**JAK1 is a genuine check on the scaling.** Three pilots at different depths, using
different cells and different candidate sets, agree to within a factor of two: 10, 12 and
19 cells per group for a ceiling of 0.7. Nothing forced that agreement.

## Caveats

**The candidate set is part of the answer.** Nearest-competitor separation depends on how
many competitors there are, so the calibration is bucketed by candidate count and a
prescription is only valid for the set size it was computed at. The GATA1 rows illustrate
the cost of ignoring this: the eligible variant count falls from 247 to 209 to 107 as depth
rises, so the three pilots are not measuring the same benchmark, and only the deepest
candidate set clears the gate.

**The dialled configurations are not biology.** Their directions come from real cell means
and their noise is real, but their rank and amplitude spread are imposed. They establish
which measurement statistic governs discrimination; they do not claim that any real dataset
has a particular geometry.

**The calibration is monotone by construction.** A binned median with a running maximum
assumes only that the ceiling rises with resolution, which is what the sweep supports. It
is not a functional form and should not be read as one.

**The nearest-competitor estimator is biased upward near the floor.** In a controlled check
a true separation of 0.029 was estimated at 0.056. It ranks configurations reliably and is
not accurate as a point estimate where the separation is far below the noise, which is the
regime the refusals fall in; the refusals do not depend on that accuracy, only on the sign
of the lower bound.
