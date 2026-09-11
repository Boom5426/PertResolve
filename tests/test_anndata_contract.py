"""Regression tests for the explicit AnnData preprocessing contract."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pertresolve.resolution import preprocess_anndata, profiles_from_anndata


def fake_adata(X):
    return SimpleNamespace(
        X=X,
        layers={},
        obsm={"latent": np.arange(12, dtype=np.float32).reshape(6, 2)},
        obs=pd.DataFrame({"perturbation": ["ctrl", "ctrl", "v1", "v1", "v2", "v2"]}),
    )


def test_profiles_from_anndata_records_effective_representation():
    adata = fake_adata(np.arange(18, dtype=np.float32).reshape(6, 3))
    grouped = profiles_from_anndata(
        adata,
        perturbation_key="perturbation",
        control="ctrl",
        depth=1,
        n_groups=2,
        seed=11,
        use_rep="latent",
        n_components=1,
    )

    assert grouped.profiles.shape == (2, 2, 2)
    meta = grouped.preprocessing
    assert meta["source"] == "obsm['latent']"
    assert meta["source_kind"] == "obsm"
    assert meta["reduction"] == "none"
    assert meta["requested_n_components"] == 1
    assert meta["n_components_applied"] is None
    assert meta["normalization"] == "none"
    assert meta["log1p"] is False
    assert meta["hvg_selection"] is False
    assert meta["perturbation_key"] == "perturbation"
    assert meta["sampling_seed"] == 11


def test_sparse_dense_materialisation_emits_memory_warning():
    sparse = pytest.importorskip("scipy.sparse")
    adata = fake_adata(sparse.csr_matrix(np.ones((6, 3), dtype=np.float32)))

    with pytest.warns(UserWarning, match="dense"):
        X, meta = preprocess_anndata(
            adata, n_components=None, dense_warning_gb=1e-12,
        )

    assert isinstance(X, np.ndarray)
    assert meta["sparse_input"] is True
    assert meta["dense_memory_estimate_bytes"] == 6 * 3 * 4
    assert meta["dense_memory_warning"]


def test_layer_and_representation_are_mutually_exclusive():
    adata = fake_adata(np.ones((6, 3)))
    adata.layers["counts"] = np.ones((6, 3))
    with pytest.raises(ValueError, match="mutually exclusive"):
        preprocess_anndata(adata, layer="counts", use_rep="latent")
