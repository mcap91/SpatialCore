"""Pytest configuration and fixtures for SpatialCore tests."""

import numpy as np
import pytest

import anndata as ad


@pytest.fixture
def mock_spatial_adata():
    """Create a mock spatial AnnData object for testing."""
    n_obs = 100
    n_vars = 50

    # Random expression matrix
    X = np.random.poisson(5, size=(n_obs, n_vars)).astype(np.float32)

    # Create AnnData
    adata = ad.AnnData(X)
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.var_names = [f"gene_{i}" for i in range(n_vars)]

    # Add spatial coordinates
    adata.obsm["spatial"] = np.random.rand(n_obs, 2) * 1000

    # Add some obs columns
    adata.obs["cell_type"] = np.random.choice(
        ["TypeA", "TypeB", "TypeC"], size=n_obs
    )

    return adata


@pytest.fixture
def spatial_weights_matrix(mock_spatial_adata):
    """Create a spatial weights matrix for testing."""
    from scipy.spatial import distance_matrix

    coords = mock_spatial_adata.obsm["spatial"]
    dist = distance_matrix(coords, coords)

    # k-nearest neighbors (k=6)
    k = 6
    W = np.zeros_like(dist)
    for i in range(len(dist)):
        neighbors = np.argsort(dist[i])[1:k+1]  # exclude self
        W[i, neighbors] = 1

    # Row-normalize
    W = W / W.sum(axis=1, keepdims=True)

    return W
