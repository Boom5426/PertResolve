"""Contracts the repository states in prose, checked instead of trusted.

The public repository makes a few structural guarantees that are cheap to
check automatically, so the display tree and its safety rules do not drift.

  * ``--out`` is refused if it resolves inside ``results/``. Twenty-six of the thirty-eight
    scripts taking ``--out`` never called :func:`reject_repo_results`, including
    ``scripts/run_all_splits.py``, the one the README names one line above the sentence.
  * The package version is one number. The import namespaces and distribution metadata
    must agree.
    distribution metadata said 0.2.0.
  * ``manuscript/`` contains only the three final deliverables intended for GitHub.

These tests need no external data and no compute workspace: two read source text, one
runs a script that reads only committed tables.

Runs standalone as well as under pytest:

    python tests/test_repo_contracts.py
"""

from __future__ import annotations

import ast
import csv
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from pertresolve.paths import reject_repo_results  # noqa: E402

#: Directories holding runnable analysis scripts. ``pertresolve/`` is deliberately absent:
#: the packaged CLI must keep working from a non-editable install, where the repository
#: (and therefore ``results/``) is not on disk at all and ``repo_root()`` cannot resolve.
SCRIPT_ROOTS = ("scripts", "results")


def _declares_out(tree: ast.Module) -> bool:
    """True if the module registers an ``--out`` command-line argument."""
    return any(isinstance(node, ast.Constant) and node.value == "--out"
               for node in ast.walk(tree))


def scripts_taking_out() -> list[Path]:
    """Every script under :data:`SCRIPT_ROOTS` that accepts ``--out``."""
    found = []
    for root in SCRIPT_ROOTS:
        for path in sorted((REPO / root).rglob("*.py")):
            tree = ast.parse(path.read_text(), filename=str(path))
            if _declares_out(tree):
                found.append(path)
    assert found, "no scripts with --out were discovered; the search is broken, not the repo"
    return found


def test_every_out_taking_script_guards_the_results_tree():
    """Each ``--out`` script must route the value through ``reject_repo_results``.

    Static rather than behavioural because most of these scripts cannot start without the
    compute workspace, which is not redistributed. ``test_the_guard_actually_refuses``
    below closes that gap on the one script that runs from committed tables alone.
    """
    unguarded = [p.relative_to(REPO) for p in scripts_taking_out()
                 if "reject_repo_results" not in p.read_text()]
    assert not unguarded, (
        "these scripts accept --out without guarding the committed results/ tree:\n  "
        + "\n  ".join(str(p) for p in unguarded))


def test_the_guard_actually_refuses():
    """The static check above means nothing if the guard itself stopped refusing."""
    for candidate in (REPO / "results", REPO / "results" / "canonical" / "x.csv"):
        try:
            reject_repo_results(candidate)
        except SystemExit as exc:
            assert "may not point inside" in str(exc)
        else:                                    # pragma: no cover - only on a regression
            raise AssertionError(f"{candidate} was accepted")
    outside = reject_repo_results(REPO / "build" / "scratch")
    assert outside == (REPO / "build" / "scratch").resolve()


def test_a_real_script_refuses_a_results_path():
    """End to end, on the one script that needs no external data."""
    script = REPO / "results" / "pilot_validation" / "pilot_validate.py"
    proc = subprocess.run(
        [sys.executable, str(script), "--out", str(REPO / "results" / "pilot_validation")],
        capture_output=True, text=True, cwd=REPO)
    assert proc.returncode != 0, "writing into results/ succeeded"
    assert "may not point inside" in proc.stdout + proc.stderr


def test_version_is_declared_once():
    """The importable package must equal the installed PertResolve version."""
    from importlib.metadata import PackageNotFoundError, version

    import pertresolve

    try:
        installed = version("pertresolve")
    except PackageNotFoundError:                 # pragma: no cover - source-tree run
        import pytest
        pytest.skip("pertresolve is not installed; nothing to compare against")
    assert pertresolve.__version__ == installed


def test_public_manuscript_contains_only_final_deliverables():
    """Keep the GitHub-facing manuscript directory intentionally minimal."""
    names = sorted(path.name for path in (REPO / "manuscript").iterdir())
    assert names == [
        "PertResolve_SI.pdf",
        "PertResolve_manuscript.pdf",
        "supplementary_data_1_datasets.xlsx",
    ]


def test_public_benchmark_has_explicit_variant_and_reference_counts():
    """PertResolve-Bench distinguishes variant conditions from its two WT rows."""
    path = REPO / "data" / "pertresolve_bench.csv"
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 472
    assert len([row for row in rows if row["variant"] != "WT"]) == 470
    assert len([row for row in rows if row["variant"] == "WT"]) == 2
    assert "is_hotspot_leaked_OLD" not in rows[0]


def test_self_contained_resolution_demo_runs():
    """The README quick-start path must work without a hidden data workspace."""
    proc = subprocess.run(
        [sys.executable, str(REPO / "examples" / "run_resolution_demo.py")],
        capture_output=True, text=True, cwd=REPO)
    assert proc.returncode == 0, proc.stderr
    assert "matrix: 160 cells x 939 features" in proc.stdout
    assert "resolution report" in proc.stdout
    assert "detection" in proc.stdout
    assert "identification" in proc.stdout
    assert "model-ranking resolution" in proc.stdout
    assert "verdict:" not in proc.stdout


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
