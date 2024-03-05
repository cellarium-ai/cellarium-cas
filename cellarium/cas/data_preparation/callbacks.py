import typing as t

import numpy as np

if t.TYPE_CHECKING:
    import anndata

TOTAL_MRNA_UMIS_COLUMN_NAME = "total_mrna_umis"


def calculate_total_mrna_umis(adata: "anndata.AnnData") -> None:
    """
    Calculate the total mRNA UMIs (Unique Molecular Identifiers) for each observation in the AnnData object and add them
    as a new column in the ``.obs`` attribute. It is recommended to use this callback before data sanitization to
    calculate all expressed genes.

    :param adata: The annotated data matrix from which to calculate the total mRNA UMIs.

    :return: A new AnnData object with the same data as ``adata`` but with an additional field in ``.obs`` containing
        the total mRNA UMIs for each observation.

    .. note::
       The function assumes that the ``.X`` attribute of the input AnnData is a count matrix where each element i,j
       represents the count of unique transcripts (UMIs) for the ith observation and jth variable (gene). The function
       performs an inline operation on the input AnnData object, changing it in place.

    Example:
    ________
        >>> calculate_total_mrna_umis(adata=adata)
        >>> 'total_mrna_umis' in adata.obs.columns
        True.
    """
    adata.obs[TOTAL_MRNA_UMIS_COLUMN_NAME] = np.array(adata.X.sum(axis=1)).flatten()


_PRE_SANITIZE_CALLBACKS: t.List[t.Callable[["anndata.AnnData"], None]] = [calculate_total_mrna_umis]


def pre_sanitize_callback(adata: "anndata.AnnData") -> None:
    """
    Apply each necessary callback before data sanitization

    :param adata: Input :class:`anndata.AnnData` instance

    :return: A new :class:`anndata.AnnData` instance after callback modifications
    """
    for callback in _PRE_SANITIZE_CALLBACKS:
        callback(adata)
