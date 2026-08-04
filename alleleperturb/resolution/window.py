"""Whether a perturbation is separable from its control, and from its closest neighbour.

Two questions come before any model is scored, and they are not the same question.

**Detection.** Is a perturbation distinguishable from the control condition at all? This
compares the distance between two halves of the same perturbation's cells with the distance
between that perturbation and the control. If the two are comparable, the perturbation's
own sampling variation is as large as its effect and nothing can be said about it.

**Identification.** Is a perturbation distinguishable from the *other perturbations*?
Detection does not imply identification: several perturbations can each depart from the
control while remaining indistinguishable from one another, and it is identification, not
detection, that a discrimination score requires. A benchmark can pass the first and have
nothing for the second to reward.

Distances here are the energy distance between two sets of cells, which uses the whole
distribution rather than its mean. It is computed on cells, not on averaged profiles, so a
verdict does not depend on averaging having preserved the difference.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["WindowResult", "detection_window", "energy_distance"]


def energy_distance(a: np.ndarray, b: np.ndarray) -> float:
    """Energy distance between two sets of points.

    ``2 * E|A - B| - E|A - A'| - E|B - B'|``, zero when the two samples come from the same
    distribution and positive otherwise. Unlike a distance between means it responds to
    differences in spread and shape, so a perturbation that changes the variability of a
    programme without moving its average is not scored as absent.

    Args:
        a: ``(n_a, K)`` sample.
        b: ``(n_b, K)`` sample.

    Returns:
        The energy distance. Computation is quadratic in the sample sizes, which is why
        callers subsample to a fixed depth.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    ab = np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)).mean()
    aa = np.sqrt(((a[:, None, :] - a[None, :, :]) ** 2).sum(-1)).mean()
    bb = np.sqrt(((b[:, None, :] - b[None, :, :]) ** 2).sum(-1)).mean()
    return float(2 * ab - aa - bb)


@dataclass
class WindowResult:
    """Per-perturbation detection window and the rankable verdict it implies.

    Attributes:
        perturbations: names, in row order.
        d_self: median distance between two halves of the same perturbation's cells, the
            reproducibility floor.
        d_null: median distance from the perturbation to the control.
        ratio: ``d_self / d_null``. Near 1 the perturbation's own sampling variation is as
            large as its effect and the window is closed; well below 1 it is open.
        signal: ``d_null - d_self``, the margin the effect clears its own noise by.
        width: the 95% spread of ``d_self`` across splits, the margin a perturbation would
            have to clear to be called reproducibly detectable.
        rankable: ``signal > width``, per perturbation.
        n_seeds: independent splits behind each median.
        depth: cells per half.
    """

    perturbations: list[str]
    d_self: np.ndarray
    d_null: np.ndarray
    ratio: np.ndarray
    signal: np.ndarray
    width: np.ndarray
    rankable: np.ndarray
    n_seeds: int
    depth: int

    @property
    def rankable_fraction(self) -> float:
        """Share of perturbations whose effect clears their own replicate noise."""
        return float(self.rankable.mean())

    def to_frame(self):
        """A ``pandas.DataFrame``, one row per perturbation."""
        import pandas as pd
        return pd.DataFrame(dict(
            perturbation=self.perturbations, d_self=self.d_self, d_null=self.d_null,
            ratio=self.ratio, signal=self.signal, width=self.width,
            rankable=self.rankable))


def detection_window(X: np.ndarray, labels, *, control: str, depth: int = 50,
                     n_seeds: int = 8, seed: int = 0,
                     perturbations: list[str] | None = None) -> WindowResult:
    """Compare each perturbation's replicate noise with its distance from the control.

    Args:
        X: ``(n_cells, K)`` matrix.
        labels: perturbation label per cell.
        control: label of the control condition.
        depth: cells per half. The energy distance is quadratic in this, and it also sets
            what "detectable" means, since a shallower comparison detects less. It is a
            declared parameter rather than a default buried in the code.
        n_seeds: independent random halvings averaged over.
        seed: base seed.
        perturbations: restrict to these labels.

    Returns:
        A :class:`WindowResult`.

    Raises:
        ValueError: if the control label is absent, or if no perturbation has ``2 * depth``
            cells.
    """
    X = np.asarray(X, dtype=np.float64)
    labels = np.asarray(labels)
    if not np.any(labels == control):
        raise ValueError(f"control label {control!r} does not appear in labels")
    control_idx = np.where(labels == control)[0]

    if perturbations is None:
        perturbations = sorted({str(v) for v in np.unique(labels)} - {str(control)})
    usable = [p for p in perturbations if int((labels == p).sum()) >= 2 * depth]
    if not usable:
        raise ValueError(f"no perturbation has {2 * depth} cells; lower depth")

    self_draws = {p: [] for p in usable}
    null_draws = {p: [] for p in usable}
    for s in range(n_seeds):
        rng = np.random.RandomState(seed + s)
        ctrl = X[rng.choice(control_idx, min(depth, len(control_idx)), replace=False)]
        for p in usable:
            idx = np.where(labels == p)[0]
            rng.shuffle(idx)
            first, second = X[idx[:depth]], X[idx[depth:2 * depth]]
            self_draws[p].append(energy_distance(first, second))
            null_draws[p].append(energy_distance(first, ctrl))

    d_self = np.array([np.median(self_draws[p]) for p in usable])
    d_null = np.array([np.median(null_draws[p]) for p in usable])
    width = np.array([np.percentile(self_draws[p], 97.5) - np.percentile(self_draws[p], 2.5)
                      for p in usable])
    signal = d_null - d_self
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(d_null > 0, d_self / d_null, np.nan)

    return WindowResult(perturbations=usable, d_self=d_self, d_null=d_null, ratio=ratio,
                        signal=signal, width=width, rankable=signal > width,
                        n_seeds=n_seeds, depth=depth)
