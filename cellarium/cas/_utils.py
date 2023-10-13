"""
This module provides utility functions to handle matrix files and calculate optimal chunk sizes based on memory
constraints.

Constants:
-----------
RESERVE_COEFFICIENT: float
    A coefficient used to define a chunk size. This number helps to keep the chunk size below the provided limit
"""
import math
import os
import tempfile
import typing as t

import anndata

RESERVE_COEFFICIENT = 0.9


def round_down_to_nearest(x: t.Union[int, float], nearest: int) -> int:
    """
    Rounds down a number to the nearest specified number.

    :param x: The number to be rounded down.
    :param nearest: The number to which `x` should be rounded down to.
    """
    return math.floor(x / nearest) * nearest


def get_preferred_chunk_size(
    adata: "anndata.AnnData",
    chunk_size_limit_in_mbs: float = 20.0,
    round_down_to: int = 100,
) -> int:
    """
    Calculate the preferred chunk size for a given AnnData object based on its size.

    :param adata: The AnnData object for which to calculate the chunk size.
    :param chunk_size_limit_in_mbs: The upper limit of chunk size in megabytes. |br|
        `Default:` ``20.``
    :param round_down_to: The number to which the resulting chunk size should be rounded down to. |br|
        `Default:` ``200``
    :return: The preferred chunk size based on the given AnnData object's size.

    .. note::
        This function works by writing the AnnData object to a temporary file, computing the size, and then it
        calculates the chunk size by considering the specified megabytes limit and the RESERVE_COEFFICIENT.
        RESERVE_COEFFICIENT is used to be sure that chunks are less than `chunk_size_limit_in_mbs`
        Finally, it rounds the calculated chunk size to the nearest specified number.

    Example usage:

    .. code-block:: python

        adata = anndata.AnnData(np.random.randn(10000, 2000))
        chunk_size = get_preferred_chunk_size(adata=adata, chunk_size_limit_in_mbs=10, round_down_to=100)
    """
    num_cells_total = adata.shape[0]

    with tempfile.NamedTemporaryFile(suffix=".h5ad") as temp_file:
        adata.write(temp_file.name, compression="gzip")
        mb_size_adata = os.path.getsize(temp_file.name) / (1024 * 1024)

    mb_size_1_cell = mb_size_adata / num_cells_total
    chunk_size = chunk_size_limit_in_mbs / mb_size_1_cell
    chunk_size *= RESERVE_COEFFICIENT
    return round_down_to_nearest(chunk_size, nearest=round_down_to)
