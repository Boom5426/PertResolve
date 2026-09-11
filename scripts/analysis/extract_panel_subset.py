"""Stream one panel dataset's selected cells out of a resource too large to load.

Frozen protocol: docs/PREREG_RESOLUTION_PANEL_v2.md. This step selects cells and units; it
computes nothing. Selection follows D1 (unit is the target), D3 (stratum), the 200-cell
inclusion rule that is the criterion's own requirement, and the section 13c candidate-pool
cap. It never invents a rule of its own.

Two readers, because the two resources are stored differently and neither can be handed to
``anndata.read_h5ad``:

``h5ad``     an ``.h5ad`` whose X is a CSC or CSR group, read through h5py in blocks. KOLF is
             CSC over 37,567 gene columns and 7.87e9 stored values, so a row subset cannot be
             taken by indexing; the columns are streamed and filtered instead.
``xatlas``   the X-Atlas/Orion parquet mirror, whose counts live in two ragged list columns
             per cell. Read shard by shard; a pyarrow ListArray already carries the offsets a
             CSR matrix needs.
``h5ad_multi`` several ``.h5ad`` files that together hold one resource, where a stratum spans
             all of them. Tahoe-100M is 14 plate files and a cell line appears on every one,
             so the unit selection has to see all 14 before any cell is read.

Writes a directory of ``.npy`` files holding the CSR arrays, the unit label per retained
cell, and the provenance of the selection, so the run is reproducible from a small artifact
and so the consumer can memory-map the arrays instead of resident-loading them.

Usage:
  python extract_panel_subset.py --reader h5ad --source FILE --unit-key gene_target \
      --control NTC --pool-cap 1500 --out subsets/kolf
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path as _RepoPath

sys.path.insert(0, str(_RepoPath(__file__).resolve().parents[2]))
from pertresolve.paths import reject_repo_results  # noqa: E402
from pathlib import Path

import numpy as np

MIN_CELLS = 200        # n_groups * depth, the criterion's own inclusion requirement
MISSING = {"nan", "NaN", "<NA>", "NA", "None", "none", "", "-1"}


def collapse_labels(labels, rule):
    """Apply a D1 reagent-to-target collapse. Every rule is named and its groups reported."""
    if rule is None:
        return labels, {}
    if rule == "chem_formulation":
        # Tahoe stores the unit as the literal tuple "[('<drug>', <conc>, 'uM')]", measured
        # on the release. D1 asks for the target, and a chemical screen's reagent-versus-
        # target distinction is the salt, solvate or hydrate preparation, not a guide:
        # 'Belumosudil (mesylate)' and 'LY-2584702 (tosylate salt)' are preparations. A
        # TRAILING parenthetical is therefore dropped; a LEADING stereo descriptor such as
        # (S)- is part of the chemical entity and is kept, per the frozen decision that
        # stereoisomers stay separate. D2 keeps the dose, so the rebuilt unit is
        # "<entity>|<conc>" and the control becomes "dmso_tf|0.0".
        pattern = re.compile(r"^\[\('(.*)', ([^,]+), '([^']*)'\)\]$")
        out, unparsed = [], []
        for v in labels:
            m = pattern.match(str(v))
            if m is None:
                unparsed.append(str(v))
                out.append(str(v))
                continue
            name, conc, _unit = m.groups()
            entity = re.sub(r"\s*\([^()]*\)\s*$", "", name).strip().lower()
            out.append(f"{entity}|{conc.strip()}")
        if unparsed:
            raise SystemExit(
                f"{len(unparsed)} labels did not match the Tahoe unit format, for example "
                f"{unparsed[:3]}. The rule is not silently applied to labels it cannot read.")
        out = np.array(out)
        groups = {}
        for original, folded in zip(labels, out):
            groups.setdefault(folded.split("|")[0], set()).add(
                pattern.match(str(original)).group(1))
        merged = {k: sorted(v) for k, v in groups.items() if len(v) > 1}
        return out, merged
    raise SystemExit(f"unknown collapse rule {rule!r}")


def _read_codes(obs, key):
    """Return ``(codes, categories)`` without ever expanding the column to strings.

    Tahoe is 100.6 M cells across 14 files. Materialising one label column as ``<U40`` is
    about 16 GB, and the selection pass needs two columns and holds both the per-file and the
    concatenated copy, so the strings alone exhausted the node before a single count was
    read. Categorical codes are int8 to int32 and cost a few hundred megabytes for the same
    information.
    """
    import h5py

    node = obs[key]
    if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
        cats = node["categories"][:]
        cats = [c.decode() if isinstance(c, bytes) else str(c) for c in cats]
        return node["codes"][:], cats
    values = _read_col(obs, key)
    cats, codes = np.unique(values, return_inverse=True)
    return codes, [str(c) for c in cats]


def _read_col(obs, key):
    """Read one obs column as strings, resolving the categorical encoding."""
    import h5py

    node = obs[key]
    if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
        cats = node["categories"][:]
        codes = node["codes"][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cats])
        out = np.empty(codes.shape, dtype=object)
        valid = codes >= 0
        out[valid] = cats[codes[valid]]
        out[~valid] = "<NA>"
        return out.astype(str)
    arr = node[:]
    return np.array([a.decode() if isinstance(a, bytes) else str(a) for a in arr])


def select_units(labels, control, *, pool_cap, pool_seed):
    """Choose the units to score. Uniform draw, seeded, recorded."""
    names, counts = np.unique(labels, return_counts=True)
    eligible = sorted(str(n) for n, c in zip(names, counts)
                      if str(n) != control and str(n) not in MISSING and c >= MIN_CELLS)
    prov = {"n_units_eligible": len(eligible), "pool_cap": pool_cap,
            "pool_seed": pool_seed, "min_cells": MIN_CELLS}
    kept = eligible
    if pool_cap and len(eligible) > pool_cap:
        rng = np.random.default_rng(pool_seed)
        pick = np.sort(rng.choice(len(eligible), size=pool_cap, replace=False))
        kept = [eligible[i] for i in pick]
    prov["n_units_scored"] = len(kept)
    return kept, prov


def _keep_mask(labels, kept_units, control):
    wanted = set(kept_units) | {control}
    return np.fromiter((v in wanted for v in labels), dtype=bool, count=len(labels))


def read_h5ad_multi(sources, unit_key, control, stratum, pool_cap, pool_seed, block_nnz,
                    collapse, out_dir):
    """Stream several .h5ad files that together hold one resource, writing as it goes.

    Selection has to see every file before any cell is read, because a stratum spans all of
    them: a Tahoe cell line appears on all 14 plates, so a per-file selection would apply the
    200-cell rule to a fourteenth of each unit's cells.

    Two things are done in code space rather than in strings or in memory, both because the
    resource is 100.6 M cells: the selection pass carries categorical codes, and the matrix is
    written straight into memory-mapped output arrays one file at a time instead of being
    concatenated. A ``scipy.sparse.vstack`` of the finished blocks would hold the blocks and
    the result at once, which is twice a 36 GB matrix.
    """
    import h5py
    import scipy.sparse as sp

    per_file = []
    for src in sources:
        with h5py.File(src, "r") as f:
            obs = f["obs"]
            codes, cats = _read_codes(obs, unit_key)
            in_stratum = np.ones(len(codes), dtype=bool)
            for key, value in (stratum or []):
                s_codes, s_cats = _read_codes(obs, key)
                wanted = {i for i, c in enumerate(s_cats) if c == value}
                in_stratum &= np.isin(s_codes, list(wanted)) if wanted else False
            shape = tuple(int(v) for v in f["X"].attrs["shape"])
        folded, merged = collapse_labels(np.asarray(cats), collapse)
        per_file.append({"src": src, "codes": codes, "folded": folded,
                         "in_stratum": in_stratum, "shape": shape, "merged": merged})

    # Global unit vocabulary, then counts, all in code space.
    vocabulary = sorted({v for f in per_file for v in f["folded"]})
    index = {v: i for i, v in enumerate(vocabulary)}
    counts = np.zeros(len(vocabulary), dtype=np.int64)
    for f in per_file:
        f["global"] = np.array([index[v] for v in f["folded"]], dtype=np.int64)
        mapped = f["global"][f["codes"][f["in_stratum"]]]
        counts += np.bincount(mapped, minlength=len(vocabulary))

    control_index = index.get(control)
    if control_index is None:
        common = [vocabulary[i] for i in np.argsort(-counts)[:15]]
        raise SystemExit(f"control {control!r} is not a unit label. Most common: {common}")
    eligible = sorted(vocabulary[i] for i in range(len(vocabulary))
                      if i != control_index and vocabulary[i] not in MISSING
                      and counts[i] >= MIN_CELLS)
    prov = {"n_units_eligible": len(eligible), "pool_cap": pool_cap,
            "pool_seed": pool_seed, "min_cells": MIN_CELLS, "collapse": collapse,
            "collapse_groups": {k: v for f in per_file for k, v in f["merged"].items()}}
    prov["n_collapse_groups"] = len(prov["collapse_groups"])
    kept_units = eligible
    if pool_cap and len(eligible) > pool_cap:
        rng = np.random.default_rng(pool_seed)
        pick = np.sort(rng.choice(len(eligible), size=pool_cap, replace=False))
        kept_units = [eligible[i] for i in pick]
    prov["n_units_scored"] = len(kept_units)
    if len(kept_units) < 5:
        raise SystemExit(f"only {len(kept_units)} units reach {MIN_CELLS} cells")

    wanted = np.zeros(len(vocabulary), dtype=bool)
    wanted[[index[u] for u in kept_units]] = True
    wanted[control_index] = True

    # Pass one: how many rows and stored values, so the outputs can be allocated on disk.
    total_rows, total_nnz, keeps = 0, 0, []
    for f in per_file:
        keep = f["in_stratum"] & wanted[f["global"][f["codes"]]]
        keeps.append(keep)
        total_rows += int(keep.sum())
        if keep.any():
            with h5py.File(f["src"], "r") as h:
                ip = h["X"]["indptr"][:]
            rows = np.flatnonzero(keep)
            total_nnz += int((ip[rows + 1] - ip[rows]).sum())
    prov.update(n_cells_source=int(sum(len(f["codes"]) for f in per_file)),
                n_cells_kept=total_rows, kept_nnz=total_nnz, n_files=len(sources),
                encoding="csr_matrix, multi-file, streamed")

    out_dir.mkdir(parents=True, exist_ok=True)
    fmt = np.lib.format.open_memmap
    data_out = fmt(out_dir / "data.npy", mode="w+", dtype=np.float32, shape=(total_nnz,))
    idx_out = fmt(out_dir / "indices.npy", mode="w+", dtype=np.int64, shape=(total_nnz,))
    ptr_out = fmt(out_dir / "indptr.npy", mode="w+", dtype=np.int64, shape=(total_rows + 1,))
    labels_out = np.empty(total_rows, dtype=object)

    row_at, nnz_at, ptr_out[0] = 0, 0, 0
    for f, keep in zip(per_file, keeps):
        if not keep.any():
            continue
        rows = np.flatnonzero(keep)
        unit = f["global"][f["codes"][rows]]
        with h5py.File(f["src"], "r") as h:
            X = h["X"]
            ip = X["indptr"][:]
            cuts = np.flatnonzero(np.diff(rows) != 1)
            starts = np.concatenate(([0], cuts + 1))
            ends = np.concatenate((cuts + 1, [len(rows)]))
            for a, b in zip(starts, ends):
                r0, r1 = int(rows[a]), int(rows[b - 1]) + 1
                lo, hi = int(ip[r0]), int(ip[r1])
                width = hi - lo
                idx_out[nnz_at:nnz_at + width] = X["indices"][lo:hi]
                data_out[nnz_at:nnz_at + width] = X["data"][lo:hi]
                lengths = np.diff(ip[r0:r1 + 1])
                ptr_out[row_at + 1:row_at + 1 + len(lengths)] = nnz_at + np.cumsum(lengths)
                nnz_at += width
                row_at += len(lengths)
        labels_out[row_at - len(rows):row_at] = [vocabulary[u] for u in unit]
    data_out.flush(); idx_out.flush(); ptr_out.flush()
    np.save(out_dir / "shape.npy", np.asarray([total_rows, per_file[0]["shape"][1]]))
    np.save(out_dir / "labels.npy", np.asarray(labels_out, dtype=str))
    prov["n_control_cells"] = int((np.asarray(labels_out, dtype=str) == control).sum())
    (out_dir / "provenance.json").write_text(json.dumps(prov, indent=2))
    return None, None, prov


def read_h5ad(source, unit_key, control, stratum, pool_cap, pool_seed, block_nnz):
    """Stream a CSC or CSR .h5ad and return (csr, labels, provenance)."""
    import h5py
    import scipy.sparse as sp

    with h5py.File(source, "r") as f:
        obs = f["obs"]
        labels = _read_col(obs, unit_key)
        n_cells = len(labels)
        in_stratum = np.ones(n_cells, dtype=bool)
        for key, value in (stratum or []):
            in_stratum &= _read_col(obs, key) == value

        kept_units, prov = select_units(labels[in_stratum], control,
                                        pool_cap=pool_cap, pool_seed=pool_seed)
        keep = _keep_mask(labels, kept_units, control) & in_stratum
        prov.update(n_cells_source=int(n_cells), n_cells_kept=int(keep.sum()),
                    n_control_cells=int((labels[keep] == control).sum()))
        if prov["n_units_scored"] < 5:
            raise SystemExit(f"only {prov['n_units_scored']} units reach {MIN_CELLS} cells")

        X = f["X"]
        enc = X.attrs.get("encoding-type", "")
        shape = tuple(int(v) for v in X.attrs["shape"])
        indptr, indices, data = X["indptr"], X["indices"], X["data"]
        new_row = np.full(shape[0], -1, dtype=np.int64)
        new_row[np.flatnonzero(keep)] = np.arange(int(keep.sum()))
        prov["encoding"] = enc
        prov["source_nnz"] = int(data.shape[0])

        ip = indptr[:]
        rows_acc, cols_acc, vals_acc = [], [], []
        # One sequential pass. For CSC the outer axis is genes and the stored index is the
        # cell, which is why a row subset cannot be indexed and must be filtered instead.
        outer_n = len(ip) - 1
        start_outer = 0
        while start_outer < outer_n:
            end_outer = int(np.searchsorted(ip, ip[start_outer] + block_nnz, "right")) - 1
            end_outer = max(end_outer, start_outer + 1)
            end_outer = min(end_outer, outer_n)
            lo, hi = int(ip[start_outer]), int(ip[end_outer])
            if hi > lo:
                idx = indices[lo:hi]
                val = data[lo:hi]
                outer = np.searchsorted(ip, np.arange(lo, hi), "right") - 1
                if enc == "csc_matrix":
                    cell, gene = idx, outer
                else:
                    cell, gene = outer, idx
                mapped = new_row[cell]
                sel = mapped >= 0
                if sel.any():
                    rows_acc.append(mapped[sel].astype(np.int32))
                    cols_acc.append(np.asarray(gene)[sel].astype(np.int32))
                    vals_acc.append(val[sel].astype(np.float32))
            start_outer = end_outer

        rows = np.concatenate(rows_acc) if rows_acc else np.zeros(0, np.int32)
        cols = np.concatenate(cols_acc) if cols_acc else np.zeros(0, np.int32)
        vals = np.concatenate(vals_acc) if vals_acc else np.zeros(0, np.float32)
        del rows_acc, cols_acc, vals_acc
        csr = sp.coo_matrix((vals, (rows, cols)),
                            shape=(int(keep.sum()), shape[1])).tocsr()
        prov["kept_nnz"] = int(csr.nnz)
        return csr, labels[keep], prov


def read_xatlas(shards, unit_key, control, pool_cap, pool_seed, n_features):
    """Stream the X-Atlas parquet mirror. Two passes: labels, then the kept rows."""
    import pyarrow.parquet as pq
    import scipy.sparse as sp

    labels_per_shard = []
    for path in shards:
        labels_per_shard.append(
            pq.read_table(path, columns=[unit_key])[unit_key].to_numpy(zero_copy_only=False)
            .astype(str))
    labels = np.concatenate(labels_per_shard)
    kept_units, prov = select_units(labels, control, pool_cap=pool_cap, pool_seed=pool_seed)
    wanted = set(kept_units) | {control}
    prov.update(n_cells_source=int(len(labels)), n_shards=len(shards))
    if prov["n_units_scored"] < 5:
        raise SystemExit(f"only {prov['n_units_scored']} units reach {MIN_CELLS} cells")

    len_acc, idx_acc, val_acc, lab_acc = [], [], [], []
    for path, lab in zip(shards, labels_per_shard):
        sel = np.fromiter((v in wanted for v in lab), dtype=bool, count=len(lab))
        if not sel.any():
            continue
        tab = pq.read_table(path, columns=["gene_token_id", "gene_expression"])
        gt = tab["gene_token_id"].combine_chunks()
        ge = tab["gene_expression"].combine_chunks()
        offs = np.asarray(gt.offsets, dtype=np.int64)
        gvals = np.asarray(gt.values)
        evals = np.asarray(ge.values)

        # Vectorised ragged gather: build the flat positions of every retained row's run in
        # one shot rather than slicing per row. At 600,000 retained rows over 223 shards the
        # per-row loop is the whole cost, and this removes it.
        rows = np.flatnonzero(sel)
        starts = offs[rows]
        lens = offs[rows + 1] - starts
        total = int(lens.sum())
        if total:
            out_start = np.concatenate(([0], np.cumsum(lens)[:-1]))
            flat = np.repeat(starts - out_start, lens) + np.arange(total, dtype=np.int64)
            idx_acc.append(gvals[flat].astype(np.int32))
            val_acc.append(evals[flat].astype(np.float32))
        len_acc.append(lens)
        lab_acc.append(lab[sel])
        del tab, gt, ge, gvals, evals

    data = np.concatenate(val_acc) if val_acc else np.zeros(0, np.float32)
    indices = np.concatenate(idx_acc) if idx_acc else np.zeros(0, np.int32)
    all_lens = np.concatenate(len_acc) if len_acc else np.zeros(0, np.int64)
    indptr = np.concatenate(([0], np.cumsum(all_lens))).astype(np.int64)
    kept_labels = np.concatenate(lab_acc)
    csr = sp.csr_matrix((data, indices, indptr), shape=(len(kept_labels), n_features))
    prov.update(n_cells_kept=int(csr.shape[0]), kept_nnz=int(csr.nnz),
                n_control_cells=int((kept_labels == control).sum()), encoding="parquet_ragged")
    return csr, kept_labels, prov


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reader", choices=("h5ad", "h5ad_multi", "xatlas"), required=True)
    ap.add_argument("--source", required=True,
                    help="file for h5ad, glob for h5ad_multi and for xatlas shards")
    ap.add_argument("--collapse", default=None, choices=("chem_formulation",),
                    help="D1 reagent-to-target rule; the groups it forms are recorded")
    ap.add_argument("--unit-key", required=True)
    ap.add_argument("--control", required=True)
    ap.add_argument("--stratum-key", action="append", default=None)
    ap.add_argument("--stratum", action="append", default=None)
    ap.add_argument("--pool-cap", type=int, default=1500)
    ap.add_argument("--pool-seed", type=int, default=0)
    ap.add_argument("--block-nnz", type=int, default=100_000_000)
    ap.add_argument("--n-features", type=int, default=38606, help="xatlas gene-token space")
    ap.add_argument("--out", type=Path, required=True,
                    help="destination .npz; refused if it resolves inside the repository's "
                         "results/ tree")
    a = ap.parse_args(argv)
    # Same guard as every other --out script: the committed results/ tree is promoted into
    # deliberately, never written by a re-run.
    a.out = reject_repo_results(a.out)

    t0 = time.time()
    strat = list(zip(a.stratum_key or [], a.stratum or []))
    if a.reader == "h5ad":
        csr, labels, prov = read_h5ad(a.source, a.unit_key, a.control, strat,
                                      a.pool_cap, a.pool_seed, a.block_nnz)
    elif a.reader == "h5ad_multi":
        import glob
        files = sorted(glob.glob(a.source))
        if not files:
            raise SystemExit(f"no files matched {a.source!r}")
        csr, labels, prov = read_h5ad_multi(files, a.unit_key, a.control, strat,
                                            a.pool_cap, a.pool_seed, a.block_nnz,
                                            a.collapse, a.out)
        if csr is None:                       # already written, streamed
            print(json.dumps(prov, indent=2))
            return 0
    else:
        import glob
        shards = sorted(glob.glob(a.source))
        if not shards:
            raise SystemExit(f"no shards matched {a.source!r}")
        csr, labels, prov = read_xatlas(shards, a.unit_key, a.control,
                                        a.pool_cap, a.pool_seed, a.n_features)

    prov.update(source=a.source, unit_key=a.unit_key, control=a.control,
                stratum=list(zip(a.stratum_key or [], a.stratum or [])),
                seconds=round(time.time() - t0, 1))
    # A directory of plain .npy files, not an .npz. The largest stratum's three arrays are
    # 58.7 GB, and an archive has to be read whole before anything can be used, which peaks
    # near 100 GB. Separate .npy files can be memory-mapped and paged in one row block at a
    # time by the consumer.
    a.out.mkdir(parents=True, exist_ok=True)
    np.save(a.out / "data.npy", csr.data)
    np.save(a.out / "indices.npy", csr.indices)
    np.save(a.out / "indptr.npy", csr.indptr)
    np.save(a.out / "shape.npy", np.asarray(csr.shape))
    np.save(a.out / "labels.npy", labels.astype(str))
    (a.out / "provenance.json").write_text(json.dumps(prov, indent=2))
    print(json.dumps(prov, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
