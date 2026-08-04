"""Command-line entry point: point it at a dataset, get a verdict.

    alleleperturb-resolution data.h5ad --perturbation-key perturbation --control non-targeting

Prints the report and, with ``--out``, writes the summary as JSON and the per-perturbation
table as CSV so the numbers can be checked rather than taken on trust.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .profiles import profiles_from_anndata  # noqa: F401  (imported for the error message)
from .report import resolution_report

__all__ = ["main"]


def _load(path: Path, perturbation_key: str, layer, use_rep, n_components, seed):
    """Read an ``.h5ad`` and return the matrix and the labels, reducing dimension if asked."""
    try:
        import anndata as ad
    except ImportError as exc:                                   # pragma: no cover
        raise SystemExit(
            "reading .h5ad needs anndata, which is not installed. "
            "Install it, or call alleleperturb.resolution.resolution_report directly with "
            "a matrix and a list of labels."
        ) from exc
    import numpy as np

    adata = ad.read_h5ad(path)
    if perturbation_key not in adata.obs:
        raise SystemExit(
            f"{perturbation_key!r} is not a column of obs. Available columns: "
            + ", ".join(map(str, list(adata.obs.columns)[:30])))
    labels = adata.obs[perturbation_key].astype(str).to_numpy()

    if use_rep is not None:
        return np.asarray(adata.obsm[use_rep]), labels
    raw = adata.layers[layer] if layer else adata.X
    X = np.asarray(raw.todense()) if hasattr(raw, "todense") else np.asarray(raw)
    if n_components and n_components < X.shape[1]:
        try:
            from sklearn.decomposition import PCA
        except ImportError as exc:                               # pragma: no cover
            raise SystemExit(
                "reducing dimension needs scikit-learn. Install it, pass "
                "--n-components 0 to keep the full space, or supply --use-rep."
            ) from exc
        X = PCA(n_components=n_components, random_state=seed).fit_transform(
            X.astype("float32"))
    return X, labels


def main(argv: list[str] | None = None) -> int:
    """Run the command-line interface. Returns a process exit status."""
    ap = argparse.ArgumentParser(
        prog="alleleperturb-resolution",
        description="Report whether a perturbation dataset can arbitrate a model comparison.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    ap.add_argument("h5ad", type=Path, help="AnnData file of single cells")
    ap.add_argument("--perturbation-key", required=True,
                    help="column of obs naming each cell's perturbation")
    ap.add_argument("--control", required=True,
                    help="value of that column marking control cells")
    ap.add_argument("--depth", type=int, default=50,
                    help="cells per group; four disjoint groups are cut per perturbation, "
                         "so a perturbation needs four times this many cells "
                         "(default: %(default)s)")
    ap.add_argument("--n-seeds", type=int, default=8,
                    help="independent cell splits averaged over (default: %(default)s)")
    ap.add_argument("--n-boot", type=int, default=1000,
                    help="bootstrap resamples for the recovery probabilities "
                         "(default: %(default)s)")
    ap.add_argument("--seed", type=int, default=0, help="base seed (default: %(default)s)")
    ap.add_argument("--n-components", type=int, default=50,
                    help="principal components to reduce to before measuring; 0 keeps the "
                         "full space (default: %(default)s)")
    ap.add_argument("--use-rep", default=None,
                    help="read obsm[USE_REP] instead of the expression matrix")
    ap.add_argument("--layer", default=None, help="read layers[LAYER] instead of X")
    ap.add_argument("--no-window", action="store_true",
                    help="skip the per-perturbation detection window, the slowest part; "
                         "the detectable fraction is then unavailable rather than guessed")
    ap.add_argument("--out", type=Path, default=None,
                    help="directory receiving resolution_report.json and, when the window "
                         "is computed, resolution_window.csv")
    args = ap.parse_args(argv)

    if not args.h5ad.exists():
        raise SystemExit(f"no such file: {args.h5ad}")

    X, labels = _load(args.h5ad, args.perturbation_key, args.layer, args.use_rep,
                      args.n_components or None, args.seed)

    import numpy as np
    if args.control not in set(map(str, labels)):
        counts = {}
        for value in labels:
            counts[str(value)] = counts.get(str(value), 0) + 1
        common = sorted(counts.items(), key=lambda kv: -kv[1])[:10]
        raise SystemExit(
            f"no cell has {args.perturbation_key} == {args.control!r}.\n"
            f"The ten most frequent values are:\n"
            + "\n".join(f"  {name!r}  {n:,} cells" for name, n in common)
            + "\nA control is usually the largest group, or a non-targeting or vehicle "
              "label. Pass one of these to --control.")

    try:
        report = resolution_report(X, labels, control=args.control, depth=args.depth,
                                   n_seeds=args.n_seeds, seed=args.seed, n_boot=args.n_boot,
                                   with_window=not args.no_window)
    except ValueError as exc:
        # These are the caller's problem to fix, not a crash to debug, so they are reported
        # as a message rather than as a traceback through the library.
        raise SystemExit(f"cannot measure this dataset: {exc}") from None
    print(report.summary())

    if args.out:
        args.out.mkdir(parents=True, exist_ok=True)
        payload = report.to_dict()
        payload["dataset"] = str(args.h5ad)
        payload["perturbation_key"] = args.perturbation_key
        payload["control"] = args.control
        (args.out / "resolution_report.json").write_text(json.dumps(payload, indent=2,
                                                                   default=float))
        if report.window is not None:
            report.window.to_frame().to_csv(args.out / "resolution_window.csv", index=False)
        print(f"\nwritten to {args.out}")
    return 0


if __name__ == "__main__":                                       # pragma: no cover
    sys.exit(main())
