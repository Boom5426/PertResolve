"""Location resolution for the analysis scripts.

The analysis scripts in ``scripts/`` and ``results/`` were originally run inside a
compute workspace that holds the per-gene expression arrays, the protein
embeddings and a shared evaluation harness. None of that is redistributed with
this repository (see the manuscript's Data availability), so the scripts cannot
assume a fixed location for it. This module resolves that location explicitly and
fails loudly when it cannot, instead of silently reading from a path that only
existed on the original machine.

Two distinct roots are involved and must not be confused:

``resolve_base``
    the directory of processed per-gene arrays and the shared scorer, supplied per run
    through ``--base`` or the ``ALLELEPERTURB_DATA`` environment variable. It is never
    inferred, because guessing it would turn a missing dependency into a confusing read
    error.

``repo_root``
    this repository, which is inferred from this file's location and holds the
    committed canonical result tables.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

__all__ = ["resolve_base", "repo_root", "add_harness_to_path", "require_inputs",
           "reject_repo_results"]

#: Files that identify the repository root, used to validate the inferred path.
_REPO_MARKERS = ("alleleperturb", "results", "setup.py")

#: Environment variable naming the directory of processed per-gene arrays and the shared
#: scorer. It was once named after an internal project, which meant nothing to anyone
#: outside it; the old name is still accepted and reported when used.
ENV_VAR = "ALLELEPERTURB_DATA"
LEGACY_ENV_VAR = "VCCOMPASS_BASE"

#: Environment variable naming the directory holding the external ``.h5ad``
#: perturbation atlases (Replogle, Norman, Adamson, sci-Plex, VCC). These are
#: public but too large to redistribute, and they live outside the workspace.
ATLAS_ENV_VAR = "ALLELEPERTURB_ATLAS_DIR"


def repo_root() -> Path:
    """Return this repository's root directory.

    Raises:
        RuntimeError: if the inferred directory does not look like the repository,
            which happens when this module is copied out of the package.
    """
    candidate = Path(__file__).resolve().parents[1]
    missing = [m for m in _REPO_MARKERS if not (candidate / m).exists()]
    if missing:
        raise RuntimeError(
            f"Inferred repository root {candidate} is missing {', '.join(missing)}. "
            "alleleperturb/paths.py must stay inside the repository."
        )
    return candidate


def resolve_base(cli_value: str | os.PathLike[str] | None = None,
                 *,
                 what: str = "the compute workspace",
                 env_var: str = ENV_VAR,
                 flag: str = "--base") -> Path:
    """Resolve an external directory that is not part of this repository.

    Resolution order is ``cli_value``, then ``env_var``. There is deliberately no
    fallback: these directories are not part of the repository, so any inferred
    value would be wrong and would turn a missing dependency into a confusing
    read error somewhere further down.

    Args:
        cli_value: value of the script's command-line argument, if given.
        what: short description of the directory, used in the error message.
        env_var: environment variable consulted when ``cli_value`` is absent.
        flag: the argument name to quote back to the user in errors.

    Returns:
        The resolved, existing directory.

    Raises:
        SystemExit: if neither source supplies a value, or the value does not
            point at an existing directory. Both cases name the flag and the
            environment variable so the caller knows how to proceed.
    """
    raw = cli_value or os.environ.get(env_var)
    if not raw and env_var == ENV_VAR and os.environ.get(LEGACY_ENV_VAR):
        raw = os.environ[LEGACY_ENV_VAR]
        print(f"note: {LEGACY_ENV_VAR} is the former name of {ENV_VAR} and is still "
              f"honoured; set {ENV_VAR} instead.", file=sys.stderr)
    if not raw:
        raise SystemExit(
            f"Could not locate {what}. Pass {flag} /path/to/dir or set "
            f"{env_var}. This directory is not redistributed with this "
            "repository; see the repository README for the expected layout."
        )
    base = Path(raw).expanduser().resolve()
    if not base.is_dir():
        raise SystemExit(
            f"{what.capitalize()} does not exist: {base}\n"
            f"(resolved from {flag if cli_value else env_var})"
        )
    return base


def require_inputs(*paths: Path) -> None:
    """Fail with the full resolved path of every missing input.

    Raises:
        SystemExit: listing each missing path, so a caller pointing ``--base`` at
            the wrong directory sees exactly what was expected and where.
    """
    missing = [p for p in paths if not Path(p).exists()]
    if missing:
        listing = "\n  ".join(str(p) for p in missing)
        raise SystemExit(
            f"Missing required input(s):\n  {listing}\n"
            f"Pass --base or set {ENV_VAR} to the workspace that contains them."
        )


def reject_repo_results(out: str | os.PathLike[str], *, flag: str = "--out") -> Path:
    """Refuse an output directory that lies inside this repository's ``results/``.

    The tables under ``results/`` back manuscript numbers and are promoted there
    deliberately, one at a time, after being checked. A re-run that wrote straight
    back into that tree could silently replace a cited value with the output of a
    changed script, which is the failure this guard exists to make impossible.
    Write to a scratch directory and promote explicitly instead.

    Args:
        out: the directory the caller intends to write to. It need not exist yet.
        flag: the argument name to quote back to the user in the error message.

    Returns:
        The expanded, absolute output directory.

    Raises:
        SystemExit: if the path resolves inside ``<repo>/results``.
    """
    resolved = Path(out).expanduser().resolve()
    protected = repo_root() / "results"
    if resolved == protected or protected in resolved.parents:
        raise SystemExit(
            f"{flag} may not point inside {protected}: {resolved}\n"
            "Those tables back manuscript numbers. Write to a scratch directory, "
            "compare against the committed table, then promote the file deliberately."
        )
    return resolved


def add_harness_to_path(base: Path) -> Path:
    """Make the shared evaluation harness in ``<base>/unified`` importable.

    The harness is not part of this repository. It is added to ``sys.path`` only
    after its presence is confirmed, so a missing harness reports the path it was
    looked for at rather than an opaque ImportError.

    Args:
        base: the workspace returned by :func:`resolve_base`.

    Returns:
        The ``<base>/unified`` directory that was added to ``sys.path``.

    Raises:
        SystemExit: if ``<base>/unified/harness.py`` does not exist.
    """
    unified = Path(base) / "unified"
    require_inputs(unified / "harness.py")
    if str(unified) not in sys.path:
        sys.path.insert(0, str(unified))
    return unified
