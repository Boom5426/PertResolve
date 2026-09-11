"""Kept only so that ``pip install -e .`` works with older tooling.

The package metadata lives in ``pyproject.toml``. Duplicating it here would let the two
drift, and a dependency list that disagrees with itself is worse than one that lives in a
single place, so nothing is repeated below.
"""

from setuptools import setup

setup()
