"""PertResolve: measurement-resolution evaluation for perturbation benchmarks."""

# Must equal the version in pyproject.toml. The two drifted once, 0.1.0 here against 0.2.0
# there, so `tests/test_repo_contracts.py` now compares this against the installed
# distribution metadata and fails when they disagree.
__version__ = "0.2.0"

from .bench import PertResolveBench

__all__ = ["PertResolveBench", "__version__"]
