"""Run the frozen resolution criterion for one panel dataset of PREREG_RESOLUTION_PANEL_v2.

Frozen protocol: docs/PREREG_RESOLUTION_PANEL_v2.md. Read it before reading any number this
script produces. The criterion itself is untouched: depth 50, n_seeds 8, n_boot 1000, seed 0,
n_components 50, thresholds as released. This script only decides which cells and which
candidate units are handed to it, per decisions D1 to D9 and section 13.

One dataset is one (resource, stratum) pair. Emits:
  <name>.json        the scalar report plus the frozen provenance of this run
  <name>_units.csv   per-perturbation detection and identification, so that the jointly
                     evaluable fraction is derivable without re-running

Usage:
  python resolution_panel_v2.py --h5ad FILE --name norman2019 \
      --perturbation-key perturbation --control control --out DIR
"""
from __future__ import annotations

import argparse
import gc
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import reject_repo_results  # noqa: E402

CELL_CAP_DEFAULT = 2000        # section 8a; n_seeds * n_groups * depth = 1600
SEP = "\x1f"                   # ASCII unit separator: cannot occur in a label, so a
                               # composite unit key is injective by construction
MISSING = {"nan", "NaN", "<NA>", "NA", "None", "none", "", "-1"}
MIN_CELLS = 200                # n_groups * depth; the criterion's own inclusion requirement
N_GROUPS = 4


def _labels_and_adata(path: Path, perturbation_keys: list[str], collapse: str | None):
    """Build the perturbation unit. Several keys form a composite unit, per D2.

    The parts are joined with an ASCII unit separator, which cannot appear inside a label,
    so two distinct tuples can never collide onto one key.
    """
    import anndata as ad

    adata = ad.read_h5ad(path)
    missing = [k for k in perturbation_keys if k not in adata.obs]
    if missing:
        raise SystemExit(f"{missing!r} not columns of obs. Available: "
                         + ", ".join(map(str, list(adata.obs.columns)[:40])))
    parts = [adata.obs[k].astype(str).to_numpy() for k in perturbation_keys]
    # D7: a cell whose unit key is not fully assigned is not a perturbation. Dropping it
    # here rather than letting it through matters twice over: an unassigned level would
    # otherwise be scored as a phantom perturbation, and where the missing part is a dose it
    # would also split the control across two keys and halve the reference.
    assigned = np.ones(len(parts[0]), dtype=bool)
    for part in parts:
        assigned &= ~np.isin(part, list(MISSING))
    labels = parts[0] if len(parts) == 1 else np.array(
        [SEP.join(t) for t in zip(*parts)])
    if collapse:
        labels = _collapse(labels, collapse)
    return adata, labels, assigned


def _collapse(labels: np.ndarray, rule: str) -> np.ndarray:
    """Apply a D1 reagent-to-target collapse. Every rule is named, never inferred."""
    import re

    if rule == "guide_suffix":
        # <Gene>_gBC_<digits> and <Gene>_g<digit> style reagent ids (PerturbMulti).
        return np.array([re.sub(r"_gBC_\d+$|_g\d+$", "", v) for v in labels])
    if rule == "guide_gN":
        # <Gene>g<digit> with no separator (Papalexi eccite: ATF2g1, IFNGR2g2).
        return np.array([re.sub(r"g\d+$", "", v) for v in labels])
    raise SystemExit(f"unknown collapse rule {rule!r}")


def _select(labels: np.ndarray, control: str, *, cell_cap: int | None,
            pool_cap: int | None, pool_seed: int) -> tuple[np.ndarray, dict]:
    """Choose the cells handed to the criterion. Returns a boolean mask and its provenance."""
    names, counts = np.unique(labels, return_counts=True)
    eligible = sorted(str(n) for n, c in zip(names, counts)
                      if str(n) != control and c >= MIN_CELLS)
    prov = {"n_units_total": int((names != control).sum()),
            "n_units_eligible": len(eligible)}

    if cell_cap is None and pool_cap is None:
        # No cap active: hand the criterion every cell, exactly as
        # pertresolve.resolution.cli does, and let group_profiles do the excluding and
        # report it. Pre-filtering here would fit the reduction on a different cell set and
        # would silently empty the ``excluded`` provenance the pre-registration promises.
        prov.update(n_units_scored=len(eligible), n_units_cell_capped=0, cell_cap=None,
                    n_control_cells=int((labels == control).sum()))
        return np.ones(len(labels), dtype=bool), prov

    kept_units = eligible
    if pool_cap is not None and len(eligible) > pool_cap:
        rng = np.random.default_rng(pool_seed)
        idx = rng.choice(len(eligible), size=pool_cap, replace=False)
        kept_units = sorted(eligible[i] for i in np.sort(idx))
        prov["pool_cap"] = pool_cap
        prov["pool_seed"] = pool_seed
    prov["n_units_scored"] = len(kept_units)

    keep = np.zeros(len(labels), dtype=bool)
    keep |= labels == control
    capped = 0
    for unit in kept_units:
        at = np.flatnonzero(labels == unit)
        if cell_cap is not None and len(at) > cell_cap:
            rng = np.random.default_rng(abs(hash(unit)) % (2 ** 32))
            at = np.sort(rng.choice(at, size=cell_cap, replace=False))
            capped += 1
        keep[at] = True
    prov["n_units_cell_capped"] = capped
    prov["cell_cap"] = cell_cap
    prov["n_control_cells"] = int((labels == control).sum())
    return keep, prov


class MappedCSR:
    """A CSR matrix whose arrays stay on disk, handed out one row block at a time.

    The largest stratum holds 4.89e9 stored values. scipy forces int64 indices above 2**31,
    so its three arrays are 58.7 GB, and merely loading them peaks near 100 GB because the
    reader needs a buffer beside the destination. Nothing in the reduction needs the whole
    matrix at once: it is already blocked. Memory-mapping the arrays and materialising one
    block at a time turns a resident 58.7 GB into a resident block.

    A directory holding ``data.npy``, ``indices.npy``, ``indptr.npy`` is the on-disk form.
    """

    def __init__(self, directory: Path):
        self.data = np.load(directory / "data.npy", mmap_mode="r")
        self.indices = np.load(directory / "indices.npy", mmap_mode="r")
        self.indptr = np.load(directory / "indptr.npy", mmap_mode="r")
        self.shape = tuple(int(v) for v in np.load(directory / "shape.npy"))
        self.dtype = self.data.dtype
        self.nnz = int(self.indptr[-1])

    def block(self, lo: int, hi: int):
        """Rows ``[lo, hi)`` as an in-memory CSR."""
        import scipy.sparse as sp

        a, b = int(self.indptr[lo]), int(self.indptr[hi])
        offsets = np.asarray(self.indptr[lo:hi + 1]) - a
        return sp.csr_matrix((np.asarray(self.data[a:b]), np.asarray(self.indices[a:b]),
                              offsets), shape=(hi - lo, self.shape[1]))

    def column_mean(self, bounds) -> np.ndarray:
        """Column means, accumulated over blocks so nothing large is resident."""
        total = np.zeros(self.shape[1], dtype=np.float64)
        for lo, hi in bounds:
            total += np.asarray(self.block(lo, hi).sum(axis=0)).ravel()
        return total / self.shape[0]


def _randomized_pca_sparse(X, n_components: int, seed: int, *,
                           n_oversamples: int = 10, n_iter: int = 7,
                           block_nnz: int = 200_000_000) -> np.ndarray:
    """Exact PCA scores for a sparse matrix, without ever densifying it.

    Centring a sparse matrix destroys its sparsity, so the mean is applied implicitly:
    ``(X - 1 mu^T) B = X B - 1 (mu^T B)``, and likewise for the transpose. Everything else is
    the standard randomised range finder. The result is the same subspace ``PCA`` computes;
    component signs may differ, which the criterion cannot see because it uses distances.

    This exists because ``IncrementalPCA`` is unusable at this width: on Norman (111,445
    cells, 33,694 genes) it took 17,957 s against 55 s for dense ``PCA``, a 329-fold
    slowdown, because its per-chunk SVD cost grows with the number of genes.
    """
    n, p = X.shape
    # Everything stays in the matrix's own dtype. A float64 probe against a float32 sparse
    # matrix makes scipy upcast the whole matrix, which doubles a 33 GB working set and is
    # what killed the largest stratum.
    dt = X.dtype if X.dtype == np.float32 else np.float64
    rng = np.random.RandomState(seed)
    k = min(n_components + n_oversamples, min(n, p))

    # Row blocks of a bounded number of stored values. A single ``X.T @ B`` makes scipy
    # materialise a transposed copy of the whole matrix, which on the largest stratum is a
    # second 33 GB allocation and is what the out-of-memory kill was. Accumulating over
    # blocks bounds the extra allocation to one block instead.
    ip = np.asarray(X.indptr)
    bounds, a = [], 0
    while a < n:
        b = int(np.searchsorted(ip, ip[a] + block_nnz, "right")) - 1
        b = min(max(b, a + 1), n)
        bounds.append((a, b))
        a = b

    mapped = isinstance(X, MappedCSR)
    get = X.block if mapped else (lambda lo, hi: X[lo:hi])
    mu = (X.column_mean(bounds) if mapped
          else np.asarray(X.mean(axis=0)).ravel()).astype(dt)

    def mm(B):                                    # (X - 1 mu^T) @ B
        B = np.asarray(B, dtype=dt)
        shift = mu @ B
        out = np.empty((n, B.shape[1]), dtype=dt)
        for lo, hi in bounds:
            out[lo:hi] = np.asarray(get(lo, hi) @ B) - shift[None, :]
        return out

    def rmm(B):                                   # (X - 1 mu^T)^T @ B
        B = np.asarray(B, dtype=dt)
        out = np.zeros((p, B.shape[1]), dtype=np.float64)
        column_sums = np.zeros(B.shape[1], dtype=np.float64)
        for lo, hi in bounds:
            out += np.asarray(get(lo, hi).T @ B[lo:hi], dtype=np.float64)
            column_sums += B[lo:hi].sum(axis=0, dtype=np.float64)
        out -= np.outer(mu.astype(np.float64), column_sums)
        return out.astype(dt)

    Q, _ = np.linalg.qr(mm(rng.normal(size=(p, k)).astype(dt)))
    for _ in range(n_iter):
        Z, _ = np.linalg.qr(rmm(Q))
        Q, _ = np.linalg.qr(mm(Z))
    U_small, sv, _ = np.linalg.svd(rmm(Q).T.astype(np.float64), full_matrices=False)
    return ((Q.astype(np.float64) @ U_small)[:, :n_components]
            * sv[:n_components]).astype(np.float64)


def _densify(block) -> np.ndarray:
    return np.asarray(block.todense()) if hasattr(block, "todense") else np.asarray(block)


def _reduce(adata, keep: np.ndarray, n_components: int, seed: int, *,
            dense_budget_gb: float, batch: int,
            force: str | None = None,
            rsvd_block_nnz: int = 200_000_000) -> tuple[np.ndarray, str]:
    """Reduce the retained cells to ``n_components``, in memory if that is possible.

    Two paths, and which one ran is recorded on every result:

    ``pca``   the exact path of ``pertresolve.resolution.cli``: densify everything, then
              ``PCA``. Gate 1 shows it reproduces the v1 panel bit for bit, so it is used
              whenever the dense working set fits inside ``dense_budget_gb``.
    ``ipca``  ``IncrementalPCA`` over row chunks, densifying one chunk at a time. Needed
              above that size: sci-Plex MCF7 at 344,862 cells by 30,000 genes is 83 GB dense
              and was killed by the OOM reaper on a 124 GB node, and the largest panel
              resources are an order of magnitude beyond that again. It is mean-centred like
              ``PCA`` but is not the same arithmetic, so gate 6 measures the difference on a
              dataset where both paths run.

    Returns the reduced matrix and the name of the path taken.
    """
    import scipy.sparse as sp
    from sklearn.decomposition import PCA, IncrementalPCA

    n_rows = int(keep.sum())
    n_cols = adata.X.shape[1]
    if isinstance(adata.X, MappedCSR):
        return _randomized_pca_sparse(adata.X, n_components, seed,
                                      block_nnz=rsvd_block_nnz), "rsvd"
    if not (n_components and n_components < n_cols):
        return _densify(adata.X[keep]), "none"

    dense_gb = n_rows * n_cols * 8 / 1e9
    if dense_gb <= dense_budget_gb:
        X = _densify(adata.X[keep]).astype("float32")
        return PCA(n_components=n_components, random_state=seed).fit_transform(X), "pca"

    rows = np.flatnonzero(keep)
    block = adata.X[rows]
    if sp.issparse(block) and force != "ipca":
        return _randomized_pca_sparse(block, n_components, seed,
                                      block_nnz=rsvd_block_nnz), "rsvd"

    # Dense but over budget: fall back to IncrementalPCA. Batches are sized so that no single
    # LAPACK call sees more than 2**31 elements; sci-Plex has 110,984 genes, and a 20,000-row
    # batch there is 2.22e9 elements, which overflows LAPACK's 32-bit indexing.
    step = max(min(batch, int(0.9 * 2 ** 31 // max(n_cols, 1))), n_components * 2)
    ipca = IncrementalPCA(n_components=n_components, batch_size=step)
    for start in range(0, len(rows), step):
        chunk = rows[start:start + step]
        if len(chunk) < n_components:            # a short final chunk cannot be fitted
            break
        ipca.partial_fit(_densify(adata.X[chunk]).astype("float32"))
    out = np.empty((len(rows), n_components), dtype=np.float64)
    for start in range(0, len(rows), step):
        chunk = rows[start:start + step]
        out[start:start + len(chunk)] = ipca.transform(
            _densify(adata.X[chunk]).astype("float32"))
    return out, "ipca"


def _detect_vector(report):
    """Per-perturbation detection, aligned with the identification vector."""
    window = report.window
    if window is None:
        return np.zeros(len(report.details["perturbations"]), dtype=bool)
    lookup = {n: bool(r) for n, r in zip(window.perturbations, window.rankable)}
    return np.array([lookup.get(n, False) for n in report.details["perturbations"]])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", type=Path, default=None)
    ap.add_argument("--npz", type=Path, default=None,
                    help="a subset written by extract_panel_subset.py, for resources that "
                         "cannot be handed to anndata.read_h5ad at all. Selection was already "
                         "applied at extraction and its provenance travels with the file.")
    ap.add_argument("--name", required=True, help="panel dataset name, one (resource, stratum)")
    ap.add_argument("--perturbation-key", required=True, action="append",
                    help="repeat for a composite unit, e.g. compound then dose (D2)")
    ap.add_argument("--control", required=True, action="append",
                    help="control value, repeated in the same order as --perturbation-key")
    ap.add_argument("--collapse", default=None, help="D1 reagent-to-target rule, if needed")
    ap.add_argument("--keep-unassigned", action="store_true",
                    help="do NOT apply D7. Exists only so that a gate can separate a code "
                         "change from the frozen selection rule by reproducing the exact v1 "
                         "cohort; never use it for a panel number.")
    ap.add_argument("--drop-prefix", default=None, action="append",
                    help="drop labels starting with this, for unassigned barcode codewords "
                         "that are decode failures rather than perturbations (D1/D7)")
    ap.add_argument("--stratum-key", default=None, action="append",
                    help="D3 stratum column; repeat for a composite stratum")
    ap.add_argument("--stratum", default=None, action="append",
                    help="stratum value, repeated in the same order as --stratum-key")
    ap.add_argument("--cell-cap", type=int, default=None,
                    help=f"cells per unit; the frozen value is {CELL_CAP_DEFAULT}, "
                         "omit to run uncapped (gate 1)")
    ap.add_argument("--pool-cap", type=int, default=None, help="K_compute, frozen at 1500")
    ap.add_argument("--pool-seed", type=int, default=0)
    ap.add_argument("--candidate-pool", type=int, default=None,
                    help="restrict the SCORING candidate set to K units while keeping every "
                         "cell, so the representation basis is fitted once on the full "
                         "cohort and only the pool varies. Draws are nested by "
                         "construction: K=10 is a subset of K=20 is a subset of K=40.")
    ap.add_argument("--candidate-seed", type=int, default=0,
                    help="seed of the single permutation the nested candidate pools are "
                         "taken from")
    ap.add_argument("--matched-k", default=None,
                    type=lambda v: [int(x) for x in str(v).split(",") if x.strip()],
                    help="PREREG v3 section 4c: score every eligible unit against its own "
                         "candidate set of size K, C_v(K, r) = {v} + the first K-1 entries "
                         "of a fixed permutation with v removed. Nested across K for a fixed "
                         "draw. Unlike --candidate-pool this does not shrink the scored set.")
    ap.add_argument("--draws", type=int, default=10,
                    help="PREREG v3 section 4c: candidate draws, seeds 0 to draws-1. The "
                         "panel reports the mean across draws with the 10th and 90th "
                         "percentile as the draw-to-draw spread.")
    ap.add_argument("--cohort-min-cells", type=int, default=None,
                    help="profile only units with at least this many cells. Used to hold the "
                         "cohort fixed across a depth sweep: without it a deeper run silently "
                         "drops the units that cannot supply 4*depth cells, and the sweep "
                         "would vary the cohort as well as the depth.")
    ap.add_argument("--depth", type=int, default=50)
    ap.add_argument("--n-seeds", type=int, default=8)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n-components", type=int, default=50)
    ap.add_argument("--dense-budget-gb", type=float, default=40.0,
                    help="above this dense working set, reduce incrementally (section 14 A3)")
    ap.add_argument("--rsvd-block-nnz", type=int, default=200_000_000,
                    help="stored values per row block in the sparse reduction; bounds the "
                         "peak allocation above the matrix itself")
    ap.add_argument("--reduce-batch", type=int, default=20000,
                    help="rows per chunk on the incremental path")
    ap.add_argument("--force-reduce", choices=("pca", "ipca", "rsvd"), default=None,
                    help="force a reduction path; used by gate 6 to run both on one dataset")
    ap.add_argument("--out", type=Path, required=True,
                    help="directory receiving <name>.json and <name>_units.csv; refused if it "
                         "resolves inside the repository's results/ tree")
    args = ap.parse_args(argv)
    # The committed tables under results/ back manuscript numbers and are promoted there one at
    # a time. A re-run writing straight back into that tree could replace a cited value with the
    # output of a changed script, so the destination is checked before anything is read.
    args.out = reject_repo_results(args.out)

    from pertresolve.resolution.report import resolution_report

    if len(args.control) > len(args.perturbation_key):
        raise SystemExit("--control cannot be repeated more often than --perturbation-key")
    control = SEP.join(args.control)

    t0 = time.time()
    if args.npz is not None:
        import scipy.sparse as sp

        if args.npz.is_dir():
            csr = MappedCSR(args.npz)
            labels = np.load(args.npz / "labels.npy", allow_pickle=False).astype(str)
            extract_prov = json.loads((args.npz / "provenance.json").read_text())
        else:
            # Read one array at a time and release each before the next, so the peak is one
            # copy of the matrix rather than the matrix plus the reader's buffers.
            blob = np.load(args.npz, allow_pickle=False)
            labels = blob["labels"].astype(str)
            extract_prov = json.loads(str(blob["provenance"][0]))
            shape = tuple(blob["shape"])
            data = blob["data"]; indices = blob["indices"]; indptr = blob["indptr"]
            blob.close()
            csr = sp.csr_matrix((data, indices, indptr), shape=shape, copy=False)
            del data, indices, indptr
            gc.collect()

        class _Shim:                       # only .X is ever touched downstream
            pass
        adata = _Shim()
        adata.X = csr
        assigned = np.ones(len(labels), dtype=bool)
        n_all = int(extract_prov.get("n_cells_source", len(labels)))
    else:
        if args.h5ad is None:
            raise SystemExit("give either --h5ad or --npz")
        extract_prov = None
        adata, labels, assigned = _labels_and_adata(args.h5ad, args.perturbation_key,
                                                    args.collapse)
        n_all = len(labels)

    # The control is identified by the leading keys only when fewer --control values are
    # given than --perturbation-key. D2 makes dose part of a treatment unit; a vehicle's
    # nominal dose is not a second control condition, and splitting the reference on it
    # would shrink the very pool every distance is measured against.
    n_ctrl_keys = len(args.control)
    if n_ctrl_keys < len(args.perturbation_key):
        prefix = SEP.join(args.control) + SEP
        is_control = np.char.startswith(labels.astype(str), prefix)
        labels = np.where(is_control, control, labels)
        # A vehicle cell whose nominal dose was never recorded is still a vehicle cell. Its
        # unit key is now the control label alone and is fully assigned, so it must not be
        # dropped as unassigned: doing so would silently thin the reference pool that every
        # distance in the criterion is measured against.
        assigned = assigned | is_control

    if args.stratum_key and args.npz is None:
        if len(args.stratum or []) != len(args.stratum_key):
            raise SystemExit("--stratum must be repeated once per --stratum-key")
        in_stratum = np.ones(len(labels), dtype=bool)
        for key, value in zip(args.stratum_key, args.stratum):
            if key not in adata.obs:
                raise SystemExit(f"stratum key {key!r} is not a column of obs")
            in_stratum &= adata.obs[key].astype(str).to_numpy() == value
        if not in_stratum.any():
            raise SystemExit(f"stratum {args.stratum!r} is empty")
        adata = adata[in_stratum]
        labels = labels[in_stratum]
        assigned = assigned[in_stratum]

    if args.drop_prefix:
        for pre in args.drop_prefix:
            assigned &= ~np.char.startswith(labels.astype(str), pre)

    if args.keep_unassigned:
        assigned = np.ones(len(labels), dtype=bool)
    n_unassigned = int((~assigned).sum())
    if n_unassigned:
        adata = adata[assigned]
        labels = labels[assigned]

    if control not in set(labels):
        common = sorted(zip(*np.unique(labels, return_counts=True)[::-1]), reverse=True)[:15]
        raise SystemExit(
            f"control {control.replace(SEP, ' | ')!r} does not occur in the unit labels"
            + (f" of stratum {args.stratum!r}" if args.stratum_key else "")
            + ". Most common labels: "
            + ", ".join(f"{n.replace(SEP, ' | ')} ({c})" for c, n in common))

    keep, prov = _select(labels, control, cell_cap=args.cell_cap,
                         pool_cap=args.pool_cap, pool_seed=args.pool_seed)
    if prov["n_units_scored"] < 5:
        raise SystemExit(f"{args.name}: only {prov['n_units_scored']} units reach "
                         f"{MIN_CELLS} cells; the frozen rule needs 5")

    t_load = time.time()
    budget = args.dense_budget_gb
    if args.force_reduce == "pca":
        budget = float("inf")
    elif args.force_reduce in ("ipca", "rsvd"):
        budget = -1.0
    X, reduce_path = _reduce(adata, keep, args.n_components, args.seed,
                             dense_budget_gb=budget, batch=args.reduce_batch,
                             force=args.force_reduce,
                             rsvd_block_nnz=args.rsvd_block_nnz)
    kept_labels = labels[keep]
    t_reduce = time.time()

    cohort = None
    if args.cohort_min_cells is not None:
        names, counts = np.unique(kept_labels, return_counts=True)
        cohort = sorted(str(n) for n, c in zip(names, counts)
                        if str(n) != control and c >= args.cohort_min_cells)
        if len(cohort) < 5:
            raise SystemExit(f"{args.name}: only {len(cohort)} units hold "
                             f"{args.cohort_min_cells} cells; the frozen rule needs 5")
        prov["cohort_min_cells"] = args.cohort_min_cells
        prov["n_units_in_fixed_cohort"] = len(cohort)

    candidates = cohort
    if args.candidate_pool is not None:
        names, counts = np.unique(kept_labels, return_counts=True)
        eligible = sorted(str(n) for n, c in zip(names, counts)
                          if str(n) != control and c >= MIN_CELLS)
        if len(eligible) < args.candidate_pool:
            raise SystemExit(f"{args.name}: {len(eligible)} eligible units, fewer than the "
                             f"requested candidate pool of {args.candidate_pool}")
        # One permutation, then a prefix. Every smaller K is therefore a subset of every
        # larger one, so a K sweep varies the pool and nothing else.
        order = np.random.default_rng(args.candidate_seed).permutation(len(eligible))
        candidates = sorted(eligible[i] for i in order[:args.candidate_pool])
        prov["candidate_pool"] = args.candidate_pool
        prov["candidate_seed"] = args.candidate_seed
        prov["n_units_in_cohort"] = len(eligible)

    if args.matched_k is not None:
        from pertresolve.resolution.report import resolution_profiles, score_profiles

        pset = resolution_profiles(X, kept_labels, control=control, depth=args.depth,
                                   n_seeds=args.n_seeds, seed=args.seed,
                                   perturbations=cohort)
        units = list(pset.perturbations)
        n_units = len(units)
        too_big = [k for k in args.matched_k if k > n_units]
        if too_big:
            raise SystemExit(f"{args.name}: {n_units} eligible units, fewer than the matched "
                             f"candidate sets {too_big}")
        per_draw = []
        for k_value in args.matched_k:
            for draw in range(args.draws):
                order = np.random.default_rng(draw).permutation(n_units)
                mask = np.zeros((n_units, n_units), dtype=bool)
                for v in range(n_units):
                    competitors = [j for j in order if j != v][:k_value - 1]
                    mask[v, competitors] = True
                rep = score_profiles(pset, n_boot=args.n_boot, seed=args.seed,
                                     candidate_mask=mask)
                d = rep.details
                per_draw.append({
                    "k": k_value,
                    "draw": draw,
                    "detection_fraction": rep.detection_fraction,
                    "identification_fraction": rep.identification_fraction,
                    "jointly_evaluable_fraction": float(
                        (_detect_vector(rep) & d["identified"]).mean()),
                    "split_half_reproducibility_reference":
                        rep.split_half_reproducibility_reference,
                    "rho2_nn_median": rep.rho2_nn_median,
                    "p_rho2_nn_gt1": d["p_rho2_nn_gt1"],
                    "p_correct_order": rep.model_ranking_resolution.p_correct_order,
                    "cross_seed_identification_fraction": d.get(
                        "cross_seed_identification_fraction",
                        d.get("cross_seed_identifiable_fraction")),
                    "cross_seed_rho2_nn_median": d.get("cross_seed_rho2_nn_median"),
                    "cross_seed_p_rho2_nn_gt1": d.get("cross_seed_p_rho2_nn_gt1"),
                })
        summary = {}
        for k_value in args.matched_k:
            rows = [r for r in per_draw if r["k"] == k_value]
            entry = {}
            for key in rows[0]:
                if key in ("draw", "k"):
                    continue
                vals = np.array([r[key] for r in rows if r[key] is not None], dtype=float)
                if vals.size == 0:
                    continue
                entry[key] = {"mean": float(vals.mean()),
                              "p10": float(np.percentile(vals, 10)),
                              "p90": float(np.percentile(vals, 90))}
            summary[str(k_value)] = entry
        prov["matched_k"] = args.matched_k
        prov["draws"] = args.draws
        prov["n_units_in_cohort"] = n_units
        report = score_profiles(pset, n_boot=args.n_boot, seed=args.seed)   # native, for reference
        matched = {"per_draw": per_draw, "summary": summary}
    else:
        matched = None
        report = resolution_report(X, kept_labels, control=control, depth=args.depth,
                                   n_seeds=args.n_seeds, seed=args.seed, n_boot=args.n_boot,
                                   perturbations=candidates)
    t_score = time.time()

    args.out.mkdir(parents=True, exist_ok=True)
    d = report.details
    order = {n: i for i, n in enumerate(d["perturbations"])}
    win = report.window
    rankable = {n: bool(r) for n, r in zip(win.perturbations, win.rankable)} if win else {}
    rows = ["perturbation,detection,identification,jointly_evaluable,separation,spread"]
    for name in d["perturbations"]:
        i = order[name]
        shown = name.replace(SEP, " | ")
        det = rankable.get(name)
        ident = bool(d["identified"][i])
        joint = "" if det is None else str(bool(det) and ident)
        rows.append(f"\"{shown}\",{'' if det is None else det},{ident},{joint},"
                    f"{float(d['separations'][i]):.6g},{float(d['separation_spread'][i]):.6g}")
    (args.out / f"{args.name}_units.csv").write_text("\n".join(rows) + "\n")

    det_arr = np.array([rankable.get(n, False) for n in d["perturbations"]])
    joint_fraction = (float((det_arr & d["identified"]).mean())
                      if rankable else float("nan"))

    summary = {
        "name": args.name,
        "prereg": "docs/PREREG_RESOLUTION_PANEL_v2.md",
        "source": str(args.npz or args.h5ad),
        "extract_provenance": extract_prov,
        "perturbation_key": args.perturbation_key,
        "unit_separator": "ASCII 0x1f" if len(args.perturbation_key) > 1 else None,
        "collapse": args.collapse,
        "control": args.control,
        "stratum_key": args.stratum_key, "stratum": args.stratum,
        "n_cells_file": int(n_all), "n_cells_scored": int(keep.sum()),
        "n_cells_unassigned_dropped": n_unassigned,
        "d7_applied": not args.keep_unassigned,
        "resolution_notes": list(report.notes),
        "n_perturbations": report.n_perturbations,
        "detection_fraction": report.detection_fraction,
        "identification_fraction": report.identification_fraction,
        "jointly_evaluable_fraction": joint_fraction,
        "split_half_reproducibility_reference": report.split_half_reproducibility_reference,
        "rho2": report.rho2, "rho2_nn_median": report.rho2_nn_median, "eta2": report.eta2,
        "p_correct_order": (report.model_ranking_resolution.p_correct_order
                            if report.model_ranking_resolution is not None else None),
        "n_excluded": len(report.excluded),
        "depth": args.depth, "n_seeds": args.n_seeds, "n_boot": args.n_boot,
        "seed": args.seed, "n_components": args.n_components,
        "reduce_path": reduce_path,
        "dense_working_set_gb": round(int(keep.sum()) * adata.X.shape[1] * 8 / 1e9, 2),
        "selection": prov,
        "matched_arm": matched,
        "p_rho2_nn_gt1": report.details.get("p_rho2_nn_gt1"),
        "cross_seed_identification_fraction": report.details.get(
            "cross_seed_identification_fraction",
            report.details.get("cross_seed_identifiable_fraction")),
        "cross_seed_rho2_nn_median": report.details.get("cross_seed_rho2_nn_median"),
        "cross_seed_p_rho2_nn_gt1": report.details.get("cross_seed_p_rho2_nn_gt1"),
        "seconds": {"load": round(t_load - t0, 1), "reduce": round(t_reduce - t_load, 1),
                    "score": round(t_score - t_reduce, 1)},
    }
    (args.out / f"{args.name}.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
