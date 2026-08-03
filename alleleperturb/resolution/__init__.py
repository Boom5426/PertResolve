"""Measurement-resolution estimation for perturbation benchmarks.

The question this subpackage answers is not how well a model scores, but whether the
measurement it is scored against can reward the distinction being claimed. A benchmark
whose ground truth does not reproducibly separate two perturbations cannot order two
models that differ in how well they separate them, however the models are built.

``scaling`` estimates the two quantities that decide this: the per-profile sampling noise
and the between-perturbation signal, together with the dimensionless ratio of the two.
"""

from .scaling import (
    bootstrap_delta2,
    gram_of,
    permutation_null_delta2,
    signal_noise,
    stats_from_gram,
    tie_aware_pds,
)

__all__ = [
    "bootstrap_delta2",
    "gram_of",
    "permutation_null_delta2",
    "signal_noise",
    "stats_from_gram",
    "tie_aware_pds",
]
