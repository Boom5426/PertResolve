#!/usr/bin/env python3
"""Rebuild the TP53/KRAS arrays by joining variants to processed cell barcodes.

This is intentionally a standalone rebuild script.  It never overwrites the
frozen arrays and records the exact source files and filtering decisions.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import issparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pertresolve.paths import reject_repo_results  # noqa: E402


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_one(raw: Path, gene_sets: set[str] | None) -> tuple[dict, set[str]]:
    name = raw.name.replace(".processed.matrix.mtx.gz", "")
    gene = name.rsplit("_", 1)[-1]
    raw_dir = raw.parent
    matrix_path = raw_dir / f"{name}.processed.matrix.mtx.gz"
    genes_path = raw_dir / f"{name}.processed.genes.csv.gz"
    cells_path = raw_dir / f"{name}.processed.cells.csv.gz"
    v2c_path = raw_dir / f"{name}.variants2cell.csv.gz"

    with gzip.open(matrix_path, "rt") as f:
        mat = mmread(f)
    mat = mat.toarray() if issparse(mat) else np.asarray(mat)
    mat = mat.astype(np.float32, copy=False)
    genes = pd.read_csv(genes_path, header=None).iloc[:, 0].astype(str).tolist()
    cells = pd.read_csv(cells_path, header=None).iloc[:, 0].astype(str)
    v2c = pd.read_csv(v2c_path, sep="\t", usecols=["variant", "cell"])
    v2c["cell"] = v2c["cell"].astype(str)
    if len(cells) != mat.shape[0] or len(genes) != mat.shape[1]:
        raise ValueError(f"{gene}: matrix {mat.shape}, cells {len(cells)}, genes {len(genes)}")
    if cells.duplicated().any() or v2c["cell"].duplicated().any():
        raise ValueError(f"{gene}: duplicate cell barcode in source")
    lookup = v2c.set_index("cell")["variant"]
    labels = lookup.reindex(cells).astype("string")
    if labels.isna().any():
        raise ValueError(f"{gene}: {int(labels.isna().sum())} processed cells missing from variants2cell")
    labels = labels.astype(str).to_numpy(dtype="U64")
    filtered = {k: int((labels == k).sum()) for k in ["unassigned", "multiple"]}
    keep = ~np.isin(labels, ["unassigned", "multiple"])
    X = mat[keep]
    labels = labels[keep]
    used_genes = set(genes) if gene_sets is None else gene_sets & set(genes)
    return {
        "gene": gene,
        "X_raw": X,
        "labels": labels,
        "genes": genes,
        "cells_before": int(len(cells)),
        "cells_after": int(keep.sum()),
        "filtered": filtered,
        "source": {p.name: sha256(p) for p in [matrix_path, genes_path, cells_path, v2c_path]},
    }, used_genes


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path, required=True, help="GSE161824 raw directory")
    ap.add_argument("--out", type=Path, required=True, help="new output directory")
    args = ap.parse_args()
    args.out = reject_repo_results(args.out)
    args.out.mkdir(parents=True, exist_ok=True)
    if any(args.out.iterdir()):
        raise FileExistsError(f"output directory must be empty: {args.out}")
    records: dict[str, dict] = {}
    shared: set[str] | None = None
    for gene in ["TP53", "KRAS"]:
        rec, shared = read_one(args.raw / f"GSE161824_A549_{gene}.processed.matrix.mtx.gz", shared)
        records[gene] = rec
    assert shared is not None
    shared_genes = sorted(shared)
    for rec in records.values():
        idx = {g: i for i, g in enumerate(rec["genes"])}
        X = rec["X_raw"][:, [idx[g] for g in shared_genes]]
        mu = X.mean(axis=0, dtype=np.float64)
        sd = X.std(axis=0, dtype=np.float64)
        sd[sd == 0] = 1.0
        rec["X"] = ((X - mu) / sd).astype(np.float32)
        del rec["X_raw"]
    np.savez_compressed(
        args.out / "joint_arrays.npz",
        Xtp=records["TP53"]["X"], vtp=records["TP53"]["labels"],
        Xkr=records["KRAS"]["X"], vkr=records["KRAS"]["labels"],
        shared=np.asarray(shared_genes),
    )
    manifest = {
        "method": "processed.cells barcode join; filter variant in {unassigned,multiple}; shared genes sorted; per-gene z-score",
        "shared_genes": len(shared_genes),
        "datasets": {
            g: {k: v for k, v in r.items() if k not in {"X", "labels", "genes"}}
            | {"matrix_shape_after": list(r["X"].shape), "variants": int(pd.Series(r["labels"]).nunique())}
            for g, r in records.items()
        },
    }
    (args.out / "provenance.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    print(f"saved {args.out / 'joint_arrays.npz'}")


if __name__ == "__main__":
    main()
