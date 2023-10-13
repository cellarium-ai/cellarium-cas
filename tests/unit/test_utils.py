"""
This module contains tests related to verifying the functionality in the utility module, especially
in the context of chunking matrices and analyzing their memory sizes.

Constants:
----------
DENSITY_RATIO: float
    The ratio of non-zero elements in the generated test matrix.

TEST_GENERATED_MATRIX_N_CELLS: int
    Number of rows in the generated test matrix.

TEST_GENERATED_MATRIX_N_GENES: int
    Number of columns in the generated test matrix.

GENERATED_MATRIX_CHUNK_MEMORY_SIZE_LIMIT: float
    Memory size limit for a chunk in megabytes. Used in test with synthetically generated matrix

RELATIVE_TOLERANCE_TO_COMPARE_TEST_MATRIX: float
    Tolerance used when comparing the actual and estimated matrix sizes in tests.
"""

import math
import os
import tempfile

import anndata
import numpy as np
from scipy import sparse

from cellarium.cas import _utils

np.random.seed(0)

DENSITY_RATIO = 0.1
TEST_GENERATED_MATRIX_N_CELLS = 10000
TEST_GENERATED_MATRIX_N_GENES = 36000
GENERATED_MATRIX_CHUNK_MEMORY_SIZE_LIMIT = 20  # mbs
RELATIVE_TOLERANCE_TO_COMPARE_TEST_MATRIX = 0.08


def create_test_matrix(n_cells: int, n_genes: int, density: float) -> "anndata.AnnData":
    """
    Create a random sparse matrix with the specified dimensions and density, then create a :class:`anndata.AnnData` out
    of it

    :param n_cells: Number of rows in the matrix.
    :param n_genes: Number of columns in the matrix.
    :param density: Density of the matrix.

    :return: Test :class:`anndata.AnnData` instance
    """
    max_value = 2858
    random_sparse_matrix = sparse.random(m=n_cells, n=n_genes, density=density, format="csr", dtype=np.float64)
    random_sparse_matrix = (random_sparse_matrix * max_value).astype(np.int64)
    obs = {"sample_id": [f"sample_{i}" for i in range(n_cells)]}  # sample annotations
    var = {"gene_id": [f"gene_{i}" for i in range(n_genes)]}  # gene annotations
    return anndata.AnnData(X=random_sparse_matrix, obs=obs, var=var)


def test_get_preferred_chunk_size() -> None:
    """
    Test the functionality of the :meth:`_utils.get_preferred_chunk_size` method.

    This test generates a test matrix with a specified number of cells and genes, and then determines a preferred chunk
    size based on memory constraints. The matrix is then split into chunks of the preferred size
    (which is calculated with :meth:`_utils.get_preferred_chunk_size`)  and written to
    temporary files with gzip compression. Finally, it asserts that the largest chunk written to file does not exceed
    the specified memory limit.

    Assertions:
        - The maximum size of any created chunk file should be less than or equal
          to the defined memory size limit.
    """
    adata = create_test_matrix(
        n_cells=TEST_GENERATED_MATRIX_N_CELLS, n_genes=TEST_GENERATED_MATRIX_N_GENES, density=DENSITY_RATIO
    )

    preferred_chunk_size = _utils.get_preferred_chunk_size(
        adata=adata,
        chunk_size_limit_in_mbs=GENERATED_MATRIX_CHUNK_MEMORY_SIZE_LIMIT,
    )
    number_of_chunks = math.ceil(TEST_GENERATED_MATRIX_N_CELLS / preferred_chunk_size)

    i, j = 0, preferred_chunk_size
    actual_chunk_sizes = []

    for _ in range(number_of_chunks):
        with tempfile.NamedTemporaryFile(suffix=".h5ad") as temp_file:
            adata[i:j, :].write(temp_file.name, compression="gzip")
            actual_chunk_sizes.append(os.path.getsize(temp_file.name) / (1024 * 1024))

        i = j
        j += preferred_chunk_size

    max_actual_chunk_size = max(actual_chunk_sizes)

    assert max_actual_chunk_size <= GENERATED_MATRIX_CHUNK_MEMORY_SIZE_LIMIT, (
        f"Maximum size of a chunk was {max_actual_chunk_size} MB which is larger than the "
        f"{GENERATED_MATRIX_CHUNK_MEMORY_SIZE_LIMIT} MB limit in this test case."
    )
