"""One call that says whether a perturbation dataset can arbitrate a model comparison.

The three questions below are ordered, and the order has consequences. A perturbation that
cannot be told apart from its control cannot be told apart from another perturbation; a
benchmark whose perturbations cannot be told apart from one another cannot order two models
that differ in how well they tell them apart. Reporting a model score without establishing
where a dataset stops is how a leaderboard comes to report the measurement.

    detectable      each perturbation against the control
    identifiable    each perturbation against its closest competitor
    benchmarkable   the panel as a whole against predictors of graded, known quality

The verdict names the first level that fails, and says what would change it. Where nothing
would, it says that too: no depth rescues a separation that is not there.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .profiles import GroupedProfiles, group_profiles
from .recovery import RecoveryResult, graded_predictor_scores, ordering_recovery
from .scaling import (
    gram_of,
    nearest_neighbour_rho2,
    signal_noise,
    stats_from_gram,
    summarise_separations,
    tie_aware_pds,
)
from .window import WindowResult, detection_window

__all__ = ["ResolutionReport", "resolution_report"]

#: Verdict thresholds. They are operating points on continuous quantities, not discoveries:
#: a dataset just below one is not categorically different from one just above. They are
#: named here so a reader can move them rather than having to find them in the code.
DETECTABLE_FRACTION = 0.5      # share of perturbations clearing their own replicate noise
IDENTIFIABLE_FRACTION = 0.5    # share whose nearest competitor clears the noise
BENCHMARKABLE_CEILING = 0.65   # replicate ceiling above which graded predictors order


@dataclass
class ResolutionReport:
    """What a dataset supports, and the measurements behind it.

    Attributes:
        verdict: the deepest level the dataset reaches, one of ``"not detectable"``,
            ``"detectable"``, ``"identifiable"`` or ``"benchmarkable"``.
        reason: one sentence naming what limits it.
        detectable_fraction: share of perturbations whose effect clears their replicate noise.
        identifiable_fraction: share whose nearest-competitor separation exceeds its own
            95% spread across cell splits, the same margin detection is judged on. NaN when
            fewer than three splits were run, since the spread is then not estimable.
        replicate_ceiling: discrimination a second measurement of the same perturbation
            attains, the most a benchmark scored against this measurement can reward.
        rho2: mean squared separation over pairs, relative to noise. Reported because
            Methods names it, and known not to govern discrimination on its own.
        rho2_nn_median: median nearest-competitor separation relative to noise. This is the
            quantity that orders the attainable ceiling.
        eta2: per-profile sampling noise at the depth used.
        n_perturbations: perturbations scored.
        depth: cells per group.
        excluded: perturbations left out, with the reason.
        window: the per-perturbation detection result, if computed.
        recovery: the graded-predictor result, if computed.
        details: everything else, for callers that want the raw numbers.
    """

    verdict: str
    reason: str
    detectable_fraction: float
    identifiable_fraction: float
    replicate_ceiling: float
    rho2: float
    rho2_nn_median: float
    eta2: float
    n_perturbations: int
    depth: int
    excluded: dict[str, str] = field(default_factory=dict)
    window: WindowResult | None = None
    recovery: RecoveryResult | None = None
    details: dict[str, Any] = field(default_factory=dict)

    def summary(self) -> str:
        """A short human-readable report."""
        lines = [
            f"verdict: {self.verdict}",
            f"  {self.reason}",
            f"  {self.n_perturbations} perturbations at {self.depth} cells per group"
            + (f", {len(self.excluded)} excluded" if self.excluded else ""),
            f"  detectable   {self.detectable_fraction:6.1%} of perturbations clear their replicate noise",
            f"  identifiable {self.identifiable_fraction:6.1%} are separable from their closest competitor",
            f"  ceiling      {self.replicate_ceiling:.3f} attainable by a second measurement",
            f"  rho2         {self.rho2:+.4f} over pairs, {self.rho2_nn_median:+.4f} nearest competitor",
        ]
        if self.recovery is not None:
            lines.append(f"  ordering     {self.recovery.describe()}")
        return "\n".join(lines)

    def to_dict(self) -> dict:
        """A JSON-serialisable view, without the array-valued members."""
        out = {k: v for k, v in self.__dict__.items()
               if k not in ("window", "recovery", "details")}
        out["excluded"] = dict(self.excluded)
        if self.recovery is not None:
            out["p_correct_order"] = self.recovery.p_correct_order
            out["p_correct_winner"] = self.recovery.p_correct_winner
            out["smallest_resolved_gap"] = self.recovery.smallest_resolved_gap
            out["graded_predictors_saturated"] = self.recovery.saturated
        out.update({k: v for k, v in self.details.items() if np.isscalar(v)})
        return out


def resolution_report(X: np.ndarray, labels, *, control: str, depth: int = 50,
                      n_seeds: int = 8, seed: int = 0, n_boot: int = 1000,
                      with_window: bool = True) -> ResolutionReport:
    """Measure what a dataset supports, before any model is compared on it.

    Args:
        X: ``(n_cells, K)`` expression or embedding matrix. Reduce dimension first for raw
            counts: distances in a very high-dimensional space are dominated by directions
            carrying no perturbation signal.
        labels: perturbation label per cell.
        control: label of the control condition.
        depth: cells per group. Four disjoint groups are cut per perturbation, so a
            perturbation needs ``4 * depth`` cells to be included.
        n_seeds: independent cell splits averaged over.
        seed: base seed.
        n_boot: bootstrap resamples for the ordering-recovery probabilities.
        with_window: also run the per-perturbation detection window, which is quadratic in
            ``depth`` and is the slowest part. The detectable fraction is unavailable
            without it.

    Returns:
        A :class:`ResolutionReport`.

    Raises:
        ValueError: propagated from :func:`~alleleperturb.resolution.profiles.group_profiles`
            when too few perturbations have enough cells.
    """
    X = np.asarray(X, dtype=np.float64)
    labels = np.asarray(labels)

    # Four groups: two score the ceiling and the graded predictors, two estimate the signal.
    # Keeping them disjoint is what stops one noise realisation from setting both the
    # measurement statistic and the score it is compared against.
    grouped: GroupedProfiles = group_profiles(X, labels, control=control, depth=depth,
                                              n_groups=4, seed=seed)
    n = len(grouped)

    ssw_acc = np.zeros(n)
    d2_acc = np.zeros((n, n))
    ceilings, separations, alpha_scores = [], [], []
    for s in range(n_seeds):
        g = group_profiles(X, labels, control=control, depth=depth, n_groups=4,
                           seed=seed + s, perturbations=grouped.perturbations)
        build, truth, signal = g.profiles[:, 0], g.profiles[:, 1], g.profiles[:, 2:]
        ssw, d2 = stats_from_gram(gram_of(signal), n, 2)
        ssw_acc += ssw
        d2_acc += d2
        ceilings.append(tie_aware_pds(build, truth))
        alpha_scores.append(graded_predictor_scores(build, truth))
        separations.append(signal)

    stats = signal_noise(ssw_acc / n_seeds, d2_acc / n_seeds, 2)
    # Nearest-competitor separations are formed per split against that split's own noise and
    # then averaged. Forming them from averaged profiles over-subtracts, because averaged
    # profiles carry a fraction of the noise the debiasing term describes.
    per_split_nn = [nearest_neighbour_rho2(s, stats["eta2"]) for s in separations]
    sep_by_split = np.stack([r["separations"] for r in per_split_nn])   # (n_seeds, n)
    seps = sep_by_split.mean(axis=0)
    nn = summarise_separations(seps, stats["eta2"])

    # Identification must be judged on the same terms as detection, or the ladder inverts.
    # Detection asks that a perturbation's distance from the control exceed the full 95%
    # spread of its own replicate distance. Asking only that the nearest-competitor
    # separation have a positive point estimate is a far weaker test, and produced datasets
    # reported as 6% detectable and 95% identifiable, which is impossible when detection is
    # a precondition for identification. The margin here is the same 95% spread, taken
    # across cell splits.
    if len(sep_by_split) >= 3:
        spread = (np.percentile(sep_by_split, 97.5, axis=0)
                  - np.percentile(sep_by_split, 2.5, axis=0))
        identifiable = float((seps > spread).mean())
    else:
        # With too few splits the spread is not estimable; say so rather than fall back to
        # a weaker test that would read as the same quantity.
        identifiable = float("nan")

    ceiling = float(np.mean(ceilings))
    recovery = ordering_recovery(np.mean(alpha_scores, axis=0), n_boot=n_boot, seed=seed)

    window = None
    detectable = float("nan")
    if with_window:
        window = detection_window(X, labels, control=control, depth=depth,
                                  n_seeds=min(n_seeds, 8), seed=seed,
                                  perturbations=grouped.perturbations)
        detectable = window.rankable_fraction

    verdict, reason = _verdict(detectable, identifiable, ceiling, nn["rho2_nn_median"])

    return ResolutionReport(
        verdict=verdict, reason=reason,
        detectable_fraction=detectable,
        identifiable_fraction=identifiable,
        replicate_ceiling=ceiling,
        rho2=stats["rho2"],
        rho2_nn_median=nn["rho2_nn_median"],
        eta2=stats["eta2"],
        n_perturbations=n,
        depth=depth,
        excluded=grouped.excluded,
        window=window, recovery=recovery,
        details=dict(delta2=stats["delta2"], pair2=stats["pair2"],
                     n_nonpositive=nn["n_nonpositive"],
                     control_cells_per_reference=grouped.control_cells_per_reference,
                     n_seeds=n_seeds),
    )


def _verdict(detectable: float, identifiable: float, ceiling: float,
             rho2_nn: float) -> tuple[str, str]:
    """Name the first level that fails, and say whether more cells would change it."""
    if np.isfinite(detectable) and detectable < DETECTABLE_FRACTION:
        return ("not detectable",
                f"only {detectable:.0%} of perturbations differ from the control by more "
                f"than their own replicate noise, so most cannot be scored at all")
    if not np.isfinite(identifiable):
        return ("detectable",
                "identification was not assessed: it needs at least three cell splits to "
                "estimate the margin a separation must clear")
    if identifiable < IDENTIFIABLE_FRACTION:
        deeper = ("more cells would help: the separation is positive, only under the noise"
                  if rho2_nn > 0 else
                  "more cells would not help: the nearest-competitor separation is "
                  "estimated at or below zero, and no depth recovers a separation that is "
                  "not there")
        return ("detectable",
                f"perturbations differ from the control but only {identifiable:.0%} are "
                f"separable from their closest competitor, so a discrimination score has "
                f"nothing to reward; {deeper}")
    if ceiling < BENCHMARKABLE_CEILING:
        return ("identifiable",
                f"perturbations are separable from one another, but a second measurement "
                f"of the same perturbation reaches only {ceiling:.3f}, so model scores "
                f"here are bounded by the measurement rather than by the models")
    return ("benchmarkable",
            f"a second measurement reaches {ceiling:.3f}, above the {BENCHMARKABLE_CEILING} "
            f"operating point, so differences between models can be attributed to the models")
