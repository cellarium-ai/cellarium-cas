"""
This module contains tests related to verifying the functionality in the _io module
"""
import io

import anndata
import numpy as np

from cellarium.cas._io import adata_to_bytes


def test_adata_to_bytes():
    """
    Test the :func:`cellarium.cas._io.adata_to_bytes` function for correctness.

    This test does the following:

    1. Creates a sample anndata.AnnData object.
    2. Converts it to a bytestream using the :meth:`cellarium.cas._io.adata_to_bytes` function.
    3. Writes the bytestream to an in-memory io.BytesIO stream.
    4. Loads the anndata.AnnData object back from the io.BytesIO stream.
    5. Asserts that the loaded anndata.AnnData object is equivalent to the original one.

    :raises AssertionError: If the original and loaded AnnData objects are not the same.
    """
    data = np.random.rand(100, 100)
    original_adata = anndata.AnnData(data)

    byte_stream = adata_to_bytes(original_adata)

    with io.BytesIO(byte_stream) as f:
        loaded_adata = anndata.read_h5ad(f)

    assert np.array_equal(original_adata.X, loaded_adata.X), "Original and loaded AnnData objects are not the same!"
