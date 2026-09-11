"""Turning cells into the disjoint measurement groups every resolution estimate needs.

Everything downstream asks the same question of the data: given several measurements of the
same perturbation taken on cells that do not overlap, how do they compare with measurements
of different perturbations? That requires cutting each perturbation's cells into disjoint
groups, and the cutting is where the mistakes live, so it happens once, here.

Three rules are enforced rather than left to the caller.

**Groups come from one shuffle.** Every group a perturbation contributes is cut from a
single permutation of its cells, so disjointness is a property of the construction rather
than something to check afterwards. Drawing groups from independent shuffles leaves them
overlapping by chance, and an overlap between the cells that estimate a signal and the
cells it is scored against inflates the answer.

**Every profile in one comparison subtracts the same reference.** A profile is a cell mean
minus a control mean. Two profiles sharing a control mean have its noise cancel between
them; two using different control means do not. Mixing the two conventions inside one
estimate makes the noise term and the signal term describe different quantities, which
biases their difference by twice the control noise, enough to drive a real signal below
zero. The reference used is recorded on the result.

**Perturbations that cannot supply the groups are named, not dropped.** A perturbation
without enough cells is excluded, and which ones were excluded and why is returned
alongside the profiles.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import warnings

import numpy as np

__all__ = ["GroupedProfiles", "group_profiles", "preprocess_anndata",
           "profiles_from_anndata"]

# A dense float64 array beyond this size is easy to create accidentally from a sparse
# AnnData matrix and can exceed the memory available on a laptop.  This is a warning, not a
# hard limit: callers with enough memory can still opt into the operation and the output
# records the estimate that was made.
DENSE_WARNING_GB = 2.0


@dataclass
class GroupedProfiles:
    """Disjoint per-perturbation profile groups, with the choices that produced them.

    Attributes:
        profiles: ``(n, n_groups, K)`` array; ``profiles[v, g]`` is the mean of one disjoint
            block of perturbation ``v``'s cells, minus the control reference.
        perturbations: names in the row order of ``profiles``.
        depth: cells contributing to each group.
        n_groups: groups cut per perturbation.
        control_cells_per_reference: cells behind each control mean, which can be below
            ``depth`` when the control is small.
        shared_reference: whether all groups subtracted the same control mean.
        excluded: perturbation name to reason, for everything left out.
        preprocessing: the AnnData preprocessing contract that produced the matrix, when
            this object came from :func:`profiles_from_anndata`; ``None`` for a direct matrix.
    """

    profiles: np.ndarray
    perturbations: list[str]
    depth: int
    n_groups: int
    control_cells_per_reference: int
    shared_reference: bool
    excluded: dict[str, str] = field(default_factory=dict)
    preprocessing: dict[str, object] | None = None

    def __len__(self) -> int:
        return len(self.perturbations)

    def __repr__(self) -> str:
        return (f"GroupedProfiles({len(self)} perturbations, {self.n_groups} groups of "
                f"{self.depth} cells, {self.profiles.shape[-1]} dimensions, "
                f"{len(self.excluded)} excluded)")


def group_profiles(X: np.ndarray, labels, *, control: str | None, depth: int,
                   n_groups: int = 2, seed: int = 0,
                   shared_reference: bool = True,
                   perturbations: list[str] | None = None) -> GroupedProfiles:
    """Cut each perturbation's cells into disjoint groups and average them.

    Args:
        X: ``(n_cells, K)`` expression or embedding matrix.
        labels: perturbation label per cell, length ``n_cells``.
        control: label of the control condition. ``None`` centres on the mean of all
            non-control cells instead, which is the only option when a dataset carries no
            control; the choice is recorded on the result.
        depth: cells per group. Perturbations with fewer than ``n_groups * depth`` cells
            are excluded and listed.
        n_groups: groups to cut per perturbation, at least 2 so that noise is identifiable.
        seed: seed for the per-perturbation shuffle.
        shared_reference: keep one control mean for every group. Set ``False`` only to
            reproduce a convention that gives each group its own reference, which biases
            any estimate combining within- and between-perturbation terms.
        perturbations: restrict to these labels, in this order. Defaults to every label
            except the control, sorted.

    Returns:
        A :class:`GroupedProfiles`.

    Raises:
        ValueError: if ``n_groups`` is below 2, if ``depth`` is not positive, if ``labels``
            does not match ``X``, if ``control`` is absent from ``labels``, or if fewer
            than two perturbations survive.
    """
    if n_groups < 2:
        raise ValueError("n_groups must be at least 2; a single group leaves sampling "
                         "noise and between-perturbation signal indistinguishable")
    if depth < 1:
        raise ValueError("depth must be a positive number of cells")
    X = np.asarray(X)
    labels = np.asarray(labels)
    if len(labels) != X.shape[0]:
        raise ValueError(f"labels has length {len(labels)} but X has {X.shape[0]} rows")
    if control is not None and not np.any(labels == control):
        raise ValueError(f"control label {control!r} does not appear in labels")

    rng = np.random.RandomState(seed)
    if perturbations is None:
        perturbations = sorted({str(v) for v in np.unique(labels)} - {str(control)})

    needed = n_groups * depth
    kept, excluded = [], {}
    for name in perturbations:
        count = int((labels == name).sum())
        if count == 0:
            excluded[name] = "absent from labels"
        elif count < needed:
            excluded[name] = f"{count} cells, needs {needed} for {n_groups} groups of {depth}"
        else:
            kept.append(name)
    if len(kept) < 2:
        raise ValueError(
            f"only {len(kept)} perturbations have {needed} cells; nothing to compare. "
            f"Lower depth or n_groups, or check the control label.")

    if control is not None:
        control_idx = np.where(labels == control)[0]
    else:
        control_idx = np.where(~np.isin(labels, kept + list(excluded)))[0]
        if len(control_idx) == 0:
            control_idx = np.arange(len(labels))
    rng.shuffle(control_idx)

    if shared_reference:
        per_ref = min(depth, len(control_idx))
        references = [X[control_idx[:per_ref]].mean(axis=0)] * n_groups
    else:
        per_ref = min(depth, len(control_idx) // n_groups)
        references = [X[control_idx[g * per_ref:(g + 1) * per_ref]].mean(axis=0)
                      for g in range(n_groups)]

    out = np.empty((len(kept), n_groups, X.shape[1]), dtype=np.float64)
    for row, name in enumerate(kept):
        idx = np.where(labels == name)[0]
        rng.shuffle(idx)
        for g in range(n_groups):
            block = idx[g * depth:(g + 1) * depth]
            out[row, g] = X[block].mean(axis=0) - references[g]

    return GroupedProfiles(profiles=out, perturbations=kept, depth=depth,
                           n_groups=n_groups, control_cells_per_reference=int(per_ref),
                           shared_reference=shared_reference, excluded=excluded)


def preprocess_anndata(adata, *, layer: str | None = None, use_rep: str | None = None,
                       seed: int = 0, n_components: int | None = 50,
                       dense_warning_gb: float = DENSE_WARNING_GB) -> tuple[np.ndarray, dict]:
    """Materialise the AnnData representation under an explicit preprocessing contract.

    This entry point intentionally does *not* normalize counts, apply ``log1p`` or select
    highly variable genes.  The selected ``X``/``layers[layer]``/``obsm[use_rep]`` value is
    assumed to already be the representation the caller wants to measure.  If an expression
    matrix is selected, it is densified and optionally reduced with PCA.  Supplying
    ``use_rep`` bypasses PCA because the representation is already explicit.

    The returned metadata is JSON-serialisable and records both requested and effective
    choices, the input/output shapes and the dense working-set estimate.  A visible warning
    is emitted when that estimate reaches ``dense_warning_gb``.

    Args:
        adata: an AnnData-like object with ``X``, ``layers`` and ``obsm`` attributes.
        layer: read ``adata.layers[layer]`` instead of ``adata.X``.
        use_rep: read ``adata.obsm[use_rep]``; mutually exclusive with ``layer``.
        seed: PCA random state. It is recorded even when no reduction is applied.
        n_components: PCA dimensions. ``None`` keeps all dimensions; PCA is skipped when
            this is at least the input width.
        dense_warning_gb: warning threshold for materialising the selected matrix.

    Returns:
        ``(X, preprocessing)`` where ``X`` is a dense NumPy matrix and ``preprocessing``
        records the effective contract.

    Raises:
        ValueError: if ``layer`` and ``use_rep`` are both supplied, or ``n_components`` is
            not a positive integer or ``None``.
        KeyError: if a requested layer or representation is absent.
    """
    if layer is not None and use_rep is not None:
        raise ValueError("layer and use_rep are mutually exclusive; choose one input")
    if n_components is not None and (
            isinstance(n_components, bool) or not isinstance(n_components, (int, np.integer))
            or n_components < 1):
        raise ValueError("n_components must be a positive integer or None")
    if dense_warning_gb <= 0:
        raise ValueError("dense_warning_gb must be positive")

    if use_rep is not None:
        source = f"obsm[{use_rep!r}]"
        raw = adata.obsm[use_rep]
        source_kind = "obsm"
    elif layer is not None:
        source = f"layers[{layer!r}]"
        raw = adata.layers[layer]
        source_kind = "layer"
    else:
        source = "X"
        raw = adata.X
        source_kind = "X"

    input_shape = tuple(int(v) for v in raw.shape)
    raw_dtype = np.dtype(getattr(raw, "dtype", np.float64))
    dense_bytes = int(np.prod(input_shape, dtype=np.int64)) * raw_dtype.itemsize
    dense_gb = dense_bytes / 1e9
    dense_warning = None
    if dense_gb >= dense_warning_gb:
        dense_warning = (
            f"{source} is about {dense_gb:.2f} GB when dense; AnnData preprocessing may "
            "exceed available RAM. Supply a precomputed obsm representation or reduce "
            "the input before loading if necessary.")
        warnings.warn(dense_warning, UserWarning, stacklevel=2)

    X = np.asarray(raw.todense()) if hasattr(raw, "todense") else np.asarray(raw)
    reduction = "none"
    n_components_applied = None
    if use_rep is None and n_components is not None and n_components < X.shape[1]:
        try:
            from sklearn.decomposition import PCA
        except ImportError as exc:                           # pragma: no cover
            raise ImportError(
                "reducing dimension needs scikit-learn; install it, pass "
                "n_components=None to keep the full space, or supply use_rep"
            ) from exc
        X = PCA(n_components=n_components, random_state=seed).fit_transform(
            X.astype(np.float32))
        reduction = "pca"
        n_components_applied = int(n_components)

    preprocessing = {
        "contract_version": 1,
        "source": source,
        "source_kind": source_kind,
        "layer": layer,
        "use_rep": use_rep,
        "normalization": "none",
        "log1p": False,
        "hvg_selection": False,
        "input_shape": list(input_shape),
        "output_shape": [int(v) for v in X.shape],
        "input_dtype": str(raw_dtype),
        "output_dtype": str(X.dtype),
        "sparse_input": bool(hasattr(raw, "todense")),
        "dense_memory_estimate_bytes": dense_bytes,
        "dense_memory_estimate_gb": dense_gb,
        "dense_memory_warning": dense_warning,
        "requested_n_components": (int(n_components) if n_components is not None else None),
        "n_components_applied": n_components_applied,
        "reduction": reduction,
        "pca_random_state": int(seed) if reduction == "pca" else None,
    }
    return X, preprocessing


def profiles_from_anndata(adata, *, perturbation_key: str, control: str | None,
                          depth: int, n_groups: int = 2, seed: int = 0,
                          layer: str | None = None, use_rep: str | None = None,
                          n_components: int | None = 50,
                          dense_warning_gb: float = DENSE_WARNING_GB) -> GroupedProfiles:
    """Group an AnnData's cells, optionally reducing dimension first.

    Args:
        adata: an ``AnnData``. Only ``obs``, ``X``, ``layers`` and ``obsm`` are touched.
        perturbation_key: column of ``adata.obs`` naming each cell's perturbation.
        control: value of that column marking control cells, or ``None``.
        depth: cells per group.
        n_groups: groups per perturbation.
        seed: seed for the shuffle and for the reduction.
        layer: read ``adata.layers[layer]`` instead of ``adata.X``.
        use_rep: read ``adata.obsm[use_rep]`` instead, in which case no reduction is done.
        n_components: reduce to this many principal components before grouping. Distances
            in a very high-dimensional space are dominated by the many directions that carry
            no perturbation signal, so a reduction is the default rather than an option.
            Pass ``None`` to keep the full space.
        dense_warning_gb: warn when materialising the selected AnnData input would require
            at least this many decimal gigabytes of dense memory.

    Returns:
        A :class:`GroupedProfiles`.

    Raises:
        KeyError: if ``perturbation_key`` is not a column of ``adata.obs``.
        ImportError: if a reduction is requested and scikit-learn is not installed.
    """
    if perturbation_key not in adata.obs:
        raise KeyError(f"{perturbation_key!r} is not a column of adata.obs; "
                       f"available: {list(adata.obs.columns)[:20]}")
    labels = adata.obs[perturbation_key].astype(str).to_numpy()

    X, preprocessing = preprocess_anndata(
        adata, layer=layer, use_rep=use_rep, seed=seed,
        n_components=n_components, dense_warning_gb=dense_warning_gb)
    grouped = group_profiles(X, labels, control=control, depth=depth,
                             n_groups=n_groups, seed=seed)
    grouped.preprocessing = preprocessing
    grouped.preprocessing["perturbation_key"] = perturbation_key
    grouped.preprocessing["control"] = control
    grouped.preprocessing["depth"] = int(depth)
    grouped.preprocessing["n_groups"] = int(n_groups)
    grouped.preprocessing["sampling_seed"] = int(seed)
    return grouped
