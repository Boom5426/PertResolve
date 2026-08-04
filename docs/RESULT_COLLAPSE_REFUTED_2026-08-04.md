# The one-variable collapse asserted in Methods does not hold, 2026-08-04

Methods (Analytical scaling) states that discrimination "depends on the data only through
the dimensionless group rho and the competitor count n_test". This tests that assertion
directly and refutes it. The correct controlling quantity is identified, and the estimator
needed to use it on real data is not yet finished.

Produced with `scripts/analysis/resolution_sweep.py` and the diagnostics recorded below.
Nothing in `results/canonical/` changes as a result; this is a claim that has to be
narrowed before Phase 6 writes it into the manuscript.

## What was tested

`rho2 = delta2 / (2 * eta2)` cannot be dialled by rescaling measured profiles: multiplying
a profile scales signal and noise together and the ratio does not move. The dial therefore
has to act on the signal alone. TP53 and KRAS supply the cells, since their own
between-variant signal is measured at `|rho2| < 0.005` with every interval covering zero
(`results/canonical/floor_law_v2.csv`), so they give real cells, real covariance and real
depth scaling with essentially no separation of their own. A synthetic separation of
controlled size and controlled geometry is then added, using directions estimated from a
fifth cell block that no scoring or estimation touches.

## Result: at matched rho2, the attainable ceiling ranges from 0.63 to 0.999

Thirty-six configurations at three amplitudes, with the direction set restricted to rank
`d` and given lognormal amplitude heterogeneity of standard deviation `h`. Every row below
is at `rho2` within 1% of 0.59, `n = 20`, `m = 25`:

| rank d | heterogeneity | effective rank | ceiling |
|---|---|---|---|
| full | 0 | 17.8 | 0.999 |
| 10 | 0 | 9.6 | 0.994 |
| 5 | 0 | 4.9 | 0.977 |
| 3 | 0 | 2.9 | 0.922 |
| 2 | 0 | 2.0 | 0.864 |
| 1 | 0 | 1.0 | 0.755 |
| full | 1 | 5.7 | 0.829 |
| full | 2 | 2.5 | 0.669 |
| full | 3 | 1.6 | 0.632 |

The assertion is therefore false as written. A benchmark can have the stated
signal-to-noise and be either almost perfectly resolvable or barely above chance,
depending on a property of the configuration that `rho2` does not measure.

Two candidate explanations were tested and one was rejected:

* **A gene-shared component is not the mechanism.** Adding a term identical across all
  variants, up to 100 times the signal amplitude, moved the ceiling only from 0.999 to
  0.988 and left `rho2` unchanged, as expected since a shared term contributes nothing to
  between-variant separation.
* **Effective rank alone is not the variable either.** Rank 1.98 gives 0.864 while rank
  2.50 gives 0.669, so heterogeneity acts separately from dimensionality.

## The quantity that does control it

`rho2` is built from the mean squared separation over all pairs. A discrimination score
asks whether a perturbation's own measurement is nearer than every competitor, which only
its **closest** competitor can spoil. The two come apart whenever the separations are not
concentrated: a low-rank or heterogeneous configuration can have a large mean while many
perturbations have a neighbour buried in the noise.

Over the same 36 configurations, with the geometry known exactly so no estimator stands
between the hypothesis and the answer:

| statistic | Spearman with the ceiling |
|---|---|
| `rho2`, mean over pairs | **+0.307** |
| `rho2` from the median nearest neighbour | +0.918 |
| `rho2` from the geometric-mean nearest neighbour | **+0.933** |

## What is committed, and what is not yet true

`alleleperturb.resolution.nearest_neighbour_rho2` implements the nearest-competitor
statistic with cross-fitting: the closest competitor is chosen on one measurement and its
distance evaluated on another, because choosing and measuring on the same data selects the
pair whose noise happened to be most negative. Two tests cover it, including that it
separates a rank-1 configuration from a spread one by more than tenfold at identical mean
separation.

**It has not yet been shown to work on the real datasets.** A first attempt fed it profiles
averaged over cell splits while passing the single-split noise level, so the debiasing
over-subtracted, every separation went negative and the statistic floored at zero. The
statistic has to be formed per split and averaged, not formed from averaged profiles. Until
that is done and checked, no claim about where the four real datasets fall on the corrected
axis is supported.

Known limitation to state wherever the statistic is used: debiasing individual nearest
separations can send them negative when the true separation is far below the noise, and at
the small end the estimator is biased upward (in a controlled check, a true 0.029 was
estimated at 0.056). It ranks configurations reliably; it is not accurate as a point
estimate near the floor.

## Consequence for the manuscript

The Analytical scaling paragraph in Methods must be narrowed. What the data support is
that discrimination depends on the separation-to-noise ratio **of the nearest competitor**,
together with the competitor count, and not on the mean over pairs. The Spearman of 0.987
between `rho2` and the ceiling over the sixteen real gene-by-depth points
(`docs/RESULT_RESOLUTION_SCALING_2026-08-03.md`) is consistent with that: it is an
empirical relationship across four datasets whose geometry happens to co-vary with their
mean separation, and it should be reported as such rather than as a law.

This also bears on the design calculator. A calculator inverting a law in `rho2` alone
would give confidently wrong depths for a dataset whose configuration geometry differs
from these four, which is exactly the case it would be used for.
