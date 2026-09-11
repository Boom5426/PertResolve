"""Whether a benchmark can order predictors whose true ranking is already known.

A model score is only interpretable if the benchmark could have detected a better model.
That is a property of the measurement, not of any model, and it can be established without
running one: build predictors whose relative quality is fixed by construction, score them
through the benchmark, and ask how often the observed scores reproduce the order that was
built in.

The predictors interpolate between a wrong answer and each perturbation's own measured
effect. Quality rises with the mixing weight by construction, so any departure from that
order in the observed scores is the benchmark failing, not the predictors.

**What the wrong answer is matters, and the obvious choice does not work.** Interpolating
toward the panel mean, which is what the earlier analysis in this repository did, produces
``(1 - w) * panel_mean + w * build``. Where the panel mean is near zero this is almost
collinear with ``build`` for every positive weight, and the discrimination score is a
cosine, which cannot see a change of scale. Every predictor above zero then scores
identically: measured on a resolvable benchmark, the seven weights returned 0.500 followed
by six values of 1.000. The gradient only appears when the panel mean is large relative to
the perturbation-specific part, which is a property of the dataset rather than of the
predictors, so the construction silently measures different things on different data.

The default here instead interpolates toward another perturbation's profile, chosen by a
fixed derangement. The wrong answer is then wrong in direction rather than in scale, which
is what the score responds to, and it carries the same magnitude and covariance as the
right one so the weight is the only thing that changes. ``toward="panel_mean"`` reproduces
the earlier construction for comparison.

    The build and evaluation profiles come from disjoint cells, so the best predictor's score is
    interpreted relative to an empirical split-half reproducibility reference rather than a
    perfect-score assumption. A predictor built from the cells it is scored against would score
    1.0 by matching its own noise and would say nothing. This controlled family is a diagnostic
    construction; it is not the paper's six-head by three-feature-space model grid.

**Ordering recovery is not monotone in benchmark quality, and that is not a defect.** At the
floor no predictor is distinguishable from another. At high resolution the best two are both
near the split-half reference and are again indistinguishable, so recovery falls. Measured on dialled
configurations it peaks around a reference level of 0.705 and declines after. Nothing here
should be inverted into a design target without checking which side of that peak the data sit
on.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .scaling import tie_aware_pds

__all__ = ["RecoveryResult", "DEFAULT_WEIGHTS", "SATURATION_TOL",
           "graded_predictor_scores", "ordering_recovery"]

#: Spread of mean scores below which the graded predictors are treated as indistinguishable.
#: Set at the resolution of a score averaged over a few dozen perturbations, well under any
#: difference a benchmark could act on.
SATURATION_TOL = 1e-3

#: Mixing weights for the graded predictors, spaced so the gaps narrow toward the top where
#: telling predictors apart is hardest. Quality is monotone in this list by construction.
DEFAULT_WEIGHTS = (0.0, 0.5, 0.7, 0.85, 0.925, 0.96, 1.0)


@dataclass
class RecoveryResult:
    """How reliably the benchmark reproduces a known predictor ordering.

    Attributes:
        weights: the mixing weights used, in increasing true quality.
        mean_score: benchmark score of each predictor, averaged over perturbations.
        p_correct_order: bootstrap probability that every adjacent pair is ordered right.
        p_correct_winner: bootstrap probability that the best predictor scores highest.
        smallest_resolved_gap: smallest adjacent quality gap the benchmark orders correctly
            with probability above 0.9, or NaN if it orders none of them. This is the
            generic API's resolution in the units of the supplied synthetic weights; it is
            a Δalpha-like construction gap, not an observed ΔPDS between paper models.
        saturated: whether every predictor scored the same, to within ``SATURATION_TOL``.
            The discrimination score is a cosine and is therefore blind to scale, so when
            perturbations are widely separated the mixed prediction stays almost collinear
            with the unmixed one and every weight scores alike. The ordering probabilities
            are then zero for a reason that says nothing about the benchmark, and reading
            them as a failure would be wrong.
        tied_adjacent_pairs: indices ``k`` whose predictors ``k`` and ``k + 1`` score the
            same to within ``SATURATION_TOL``. This catches the case ``saturated`` misses
            and which is the common one at high resolution: only the *top* predictors
            collapse together while the wrong answer at weight 0 still scores far below,
            so the full range is wide and ``saturated`` stays False. ``p_correct_order``
            and ``p_correct_winner`` are then 0.00 because the comparison is a strict
            ``>`` over exact ties, which reads as the benchmark systematically choosing
            wrong when it is in fact unable to choose at all.
        n_perturbations: perturbations scored.
        n_boot: bootstrap resamples.
    """

    weights: tuple[float, ...]
    mean_score: np.ndarray
    p_correct_order: float
    p_correct_winner: float
    smallest_resolved_gap: float
    saturated: bool
    n_perturbations: int
    n_boot: int
    tied_adjacent_pairs: tuple[int, ...] = ()

    def describe(self) -> str:
        """One line stating what the benchmark can order, or why the question is empty."""
        if self.saturated:
            return (f"every graded predictor scores {self.mean_score.mean():.3f}; the "
                    f"discrimination score cannot separate them, so ordering is not assessable")
        # Reported before the probabilities, because a tie is the reason for them rather
        # than a caveat on them: a reader who sees 0.00 first has already misread it.
        tie_note = ""
        if self.tied_adjacent_pairs:
            n = len(self.tied_adjacent_pairs)
            subject = f"{n} adjacent pair of predictors scores" if n == 1 else \
                      f"{n} adjacent pairs of predictors score"
            tie_note = (f"; {subject} alike to within {SATURATION_TOL:g}, which no bootstrap "
                        f"can order, so those probabilities are limited by the tie")
        if np.isnan(self.smallest_resolved_gap):
            return (f"no adjacent pair of graded predictors is ordered reliably "
                    f"(P(full order) = {self.p_correct_order:.2f}){tie_note}")
        return (f"orders predictors differing by {self.smallest_resolved_gap:.3f} in quality "
                f"(P(full order) = {self.p_correct_order:.2f}, "
                f"P(best chosen) = {self.p_correct_winner:.2f}){tie_note}")


def graded_predictor_scores(build: np.ndarray, truth: np.ndarray,
                            weights=DEFAULT_WEIGHTS, *,
                            toward: str = "shuffled",
                            candidate_mask: np.ndarray | None = None) -> np.ndarray:
    """Score each graded predictor on each perturbation.

    Args:
        build: ``(n, K)`` profiles the predictors are built from.
        truth: ``(n, K)`` profiles they are scored against, from disjoint cells.
        weights: mixing weights, increasing.
        toward: what a weight of zero predicts. ``"shuffled"`` uses another perturbation's
            profile under a fixed derangement, so the wrong answer is wrong in direction.
            ``"panel_mean"`` reproduces the earlier construction, which a cosine score can
            only distinguish from the right answer when the panel mean is large.
        candidate_mask: optional ``(n, n)`` boolean; row ``i`` names the truths prediction
            ``i`` is ranked against. ``None`` ranks against every truth, which is the
            released behaviour. The ordering-recovery probability is a rank statistic and is
            therefore candidate-set dependent in the same way the split-half reference is, so a matched
            candidate analysis has to pass the same mask here.

    Returns:
        ``(len(weights), n)`` discrimination scores.

    Raises:
        ValueError: if ``build`` and ``truth`` disagree in shape, if fewer than two
            perturbations are supplied, or if ``toward`` is not recognised.
    """
    build = np.asarray(build, dtype=np.float64)
    truth = np.asarray(truth, dtype=np.float64)
    if build.shape != truth.shape:
        raise ValueError(f"build {build.shape} and truth {truth.shape} must match")
    n = build.shape[0]
    if n < 2:
        raise ValueError("at least two perturbations are needed to rank one against another")

    if toward == "panel_mean":
        wrong = np.broadcast_to(build.mean(axis=0), build.shape)
    elif toward == "shuffled":
        # A cyclic shift is a derangement for any n >= 2, so no perturbation is handed its
        # own profile as the wrong answer, and it needs no random state.
        wrong = np.roll(build, 1, axis=0)
    else:
        raise ValueError(f"toward must be 'shuffled' or 'panel_mean', not {toward!r}")

    out = np.empty((len(weights), n))
    for i, w in enumerate(weights):
        out[i] = _per_perturbation_pds((1.0 - w) * wrong + w * build, truth,
                                       candidate_mask=candidate_mask)
    return out


def _per_perturbation_pds(prediction: np.ndarray, truth: np.ndarray,
                          candidate_mask: np.ndarray | None = None) -> np.ndarray:
    """Tie-aware cosine discrimination score, kept per perturbation rather than averaged."""
    pn = prediction / (np.linalg.norm(prediction, axis=1, keepdims=True) + 1e-12)
    tn = truth / (np.linalg.norm(truth, axis=1, keepdims=True) + 1e-12)
    dist = 1.0 - pn @ tn.T
    n = dist.shape[1]
    scores = np.empty(dist.shape[0])
    for i in range(dist.shape[0]):
        row = dist[i]
        target = row[i]
        if candidate_mask is None:
            less = int((row < target - 1e-12).sum())
            eq = int((np.abs(row - target) <= 1e-12).sum())
            denominator = n - 1
        else:
            allowed = np.asarray(candidate_mask[i], dtype=bool).copy()
            allowed[i] = True
            sub = row[allowed]
            less = int((sub < target - 1e-12).sum())
            eq = int((np.abs(sub - target) <= 1e-12).sum())
            denominator = max(int(allowed.sum()) - 1, 1)
        scores[i] = 1.0 - (less + (eq - 1) / 2.0) / denominator
    return scores


def ordering_recovery(per_predictor: np.ndarray, weights=DEFAULT_WEIGHTS, *,
                      n_boot: int = 1000, seed: int = 0,
                      confidence: float = 0.9) -> RecoveryResult:
    """Bootstrap over perturbations to see how reliably the known order is reproduced.

    Perturbations are resampled because they are what another run of the same benchmark
    would draw differently. The question is whether this benchmark, on another panel from
    the same assay, would still order the predictors correctly.

    Args:
        per_predictor: ``(len(weights), n)`` scores from :func:`graded_predictor_scores`.
        weights: the same weights, increasing.
        n_boot: bootstrap resamples.
        seed: seed for the resampling.
        confidence: probability an adjacent pair must be ordered right for its gap to count
            as resolved.

    Returns:
        A :class:`RecoveryResult`.
    """
    per_predictor = np.asarray(per_predictor, dtype=np.float64)
    n = per_predictor.shape[1]
    rng = np.random.RandomState(seed)
    n_pairs = len(weights) - 1

    correct = winner = 0
    pair_right = np.zeros(n_pairs)
    for _ in range(n_boot):
        means = per_predictor[:, rng.randint(0, n, n)].mean(axis=1)
        ordered = np.array([means[k + 1] > means[k] for k in range(n_pairs)])
        pair_right += ordered
        correct += int(ordered.all())
        winner += int(np.argmax(means) == len(weights) - 1)
    pair_right /= n_boot

    gaps = np.array([weights[k + 1] - weights[k] for k in range(n_pairs)])
    resolved = gaps[pair_right > confidence]
    mean_score = per_predictor.mean(axis=1)
    # Ties are read off the observed means rather than the bootstrap: a pair that scores
    # identically on the data scores identically in every resample, so the bootstrap
    # cannot distinguish "always ordered wrong" from "never ordered at all".
    tied = tuple(k for k in range(n_pairs)
                 if abs(mean_score[k + 1] - mean_score[k]) < SATURATION_TOL)
    return RecoveryResult(
        weights=tuple(weights),
        mean_score=mean_score,
        p_correct_order=correct / n_boot,
        p_correct_winner=winner / n_boot,
        smallest_resolved_gap=float(resolved.min()) if len(resolved) else float("nan"),
        # ndarray.ptp() was removed in numpy 2.0; spell it out so the package works
        # on both sides of that boundary
        saturated=bool(float(mean_score.max() - mean_score.min()) < SATURATION_TOL),
        n_perturbations=n,
        n_boot=n_boot,
        tied_adjacent_pairs=tied,
    )
