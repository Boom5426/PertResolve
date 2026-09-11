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

import numpy as np

__all__ = ["GroupedProfiles", "group_profiles", "profiles_from_anndata"]


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
    """

    profiles: np.ndarray
    perturbations: list[str]
    depth: int
    n_groups: int
    control_cells_per_reference: int
    shared_reference: bool
    excluded: dict[str, str] = field(default_factory=dict)

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


def profiles_from_anndata(adata, *, perturbation_key: str, control: str | None,
                          depth: int, n_groups: int = 2, seed: int = 0,
                          layer: str | None = None, use_rep: str | None = None,
                          n_components: int | None = 50) -> GroupedProfiles:
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

    if use_rep is not None:
        X = np.asarray(adata.obsm[use_rep])
    else:
        raw = adata.layers[layer] if layer is not None else adata.X
        X = np.asarray(raw.todense()) if hasattr(raw, "todense") else np.asarray(raw)
        if n_components is not None and n_components < X.shape[1]:
            try:
                from sklearn.decomposition import PCA
            except ImportError as exc:                       # pragma: no cover
                raise ImportError(
                    "reducing dimension needs scikit-learn; install it, pass "
                    "n_components=None to keep the full space, or supply use_rep"
                ) from exc
            X = PCA(n_components=n_components, random_state=seed).fit_transform(
                X.astype(np.float32))

    return group_profiles(X, labels, control=control, depth=depth,
                          n_groups=n_groups, seed=seed)
