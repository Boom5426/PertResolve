#!/usr/bin/env python3
"""Run the self-contained PertResolve resolution demo."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DATA = REPO / "data" / "demo" / "pertresolve_resolution_demo.npz"
sys.path.insert(0, str(REPO))

from pertresolve.resolution import resolution_report  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA,
                        help="demo .npz file (default: %(default)s)")
    args = parser.parse_args()

    with np.load(args.data, allow_pickle=False) as demo:
        X = demo["X"]
        labels = demo["labels"].astype(str)
        control = str(demo["control"])
        depth = int(demo["depth"])

    report = resolution_report(
        X,
        labels,
        control=control,
        depth=depth,
        n_seeds=4,
        seed=7,
        n_boot=200,
        with_window=True,
    )
    print(f"data: {args.data}")
    print(f"matrix: {X.shape[0]:,} cells x {X.shape[1]:,} features")
    print(report.summary())


if __name__ == "__main__":
    main()
