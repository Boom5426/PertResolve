"""Measurement-resolution estimation for perturbation benchmarks.

The question here is not how well a model scores, but whether the measurement it is scored
against can reward the distinction being claimed. A benchmark whose ground truth does not
reproducibly separate two perturbations cannot order two models that differ in how well
they separate them, however the models are built, so a leaderboard below that resolution
reports the measurement as much as the methods.

Start with :func:`resolution_report`, which takes a matrix of cells and their perturbation
labels and returns a verdict, or the ``alleleperturb-resolution`` command for an ``.h5ad``.

    from alleleperturb.resolution import resolution_report
    report = resolution_report(X, labels, control="non-targeting", depth=50)
    print(report.summary())

The pieces underneath are usable on their own:

``profiles``       cut cells into the disjoint groups every estimate needs
``window``         detection, each perturbation against its control
``scaling``        sampling noise, between-perturbation signal, and their ratio
``recovery``       whether a known predictor ordering survives the benchmark
``report``         the three levels together, with a verdict

One result from calibrating these is worth knowing before reading any of them: the mean
squared separation over pairs, which earlier work in this repository used as the
measurement axis, does not govern discrimination. The nearest competitor does. At a fixed
mean separation the attainable score ranges from 0.63 to 0.999 depending on the geometry of
the configuration (``docs/RESULT_COLLAPSE_REFUTED_2026-08-04.md``), and it is
``rho2_nn_median`` rather than ``rho2`` that orders it.
"""

from .profiles import GroupedProfiles, group_profiles, profiles_from_anndata
from .recovery import (
    DEFAULT_WEIGHTS,
    SATURATION_TOL,
    RecoveryResult,
    graded_predictor_scores,
    ordering_recovery,
)
from .report import ResolutionReport, resolution_report
from .scaling import (
    bootstrap_delta2,
    gram_of,
    nearest_neighbour_rho2,
    permutation_null_delta2,
    signal_noise,
    stats_from_gram,
    summarise_separations,
    tie_aware_pds,
)
from .window import WindowResult, detection_window, energy_distance

__all__ = [
    "DEFAULT_WEIGHTS",
    "SATURATION_TOL",
    "GroupedProfiles",
    "RecoveryResult",
    "ResolutionReport",
    "WindowResult",
    "bootstrap_delta2",
    "detection_window",
    "energy_distance",
    "graded_predictor_scores",
    "gram_of",
    "group_profiles",
    "nearest_neighbour_rho2",
    "ordering_recovery",
    "permutation_null_delta2",
    "profiles_from_anndata",
    "resolution_report",
    "signal_noise",
    "stats_from_gram",
    "summarise_separations",
    "tie_aware_pds",
]
