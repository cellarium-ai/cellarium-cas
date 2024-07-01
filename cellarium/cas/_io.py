import os
import sys
import tempfile
import typing as t
import warnings
from contextlib import contextmanager

import anndata
import h5py
import numpy as np

if t.TYPE_CHECKING:
    import pathlib


def _collect_datasets(dsets: t.Dict, group: "h5py.Group"):
    for k, v in group.items():
        if isinstance(v, h5py.Dataset):
            dsets[k] = v[()]
        else:
            _collect_datasets(dsets, v)


def _read_legacy_10x_h5(filename, *, genome=None):
    """
    Read hdf5 file from Cell Ranger v2 or earlier versions.
    """
    with h5py.File(str(filename), "r") as f:
        try:
            children = list(f.keys())
            if not genome:
                if len(children) > 1:
                    raise ValueError(
                        f"'{filename}' contains more than one genome. For legacy 10x h5 "
                        "files you must specify the genome if more than one is present. "
                        f"Available genomes are: {children}"
                    )
                genome = children[0]
            elif genome not in children:
                raise ValueError(
                    f"Could not find genome '{genome}' in '{filename}'. " f"Available genomes are: {children}"
                )

            dsets = {}
            _collect_datasets(dsets, f[genome])

            # AnnData works with csr matrices
            # 10x stores the transposed data, so we do the transposition right away
            from scipy.sparse import csr_matrix

            M, N = dsets["shape"]
            data = dsets["data"]
            if dsets["data"].dtype == np.dtype("int32"):
                data = dsets["data"].view("float32")
                data[:] = dsets["data"]
            matrix = csr_matrix(
                (data, dsets["indices"], dsets["indptr"]),
                shape=(N, M),
            )
            # the csc matrix is automatically the transposed csr matrix
            # as scanpy expects it, so, no need for a further transposition
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                adata = anndata.AnnData(
                    matrix,
                    obs=dict(obs_names=dsets["barcodes"].astype(str)),
                    var=dict(
                        var_names=dsets["gene_names"].astype(str),
                        gene_ids=dsets["genes"].astype(str),
                    ),
                )

            return adata
        except KeyError:
            raise Exception("File is missing one or more required datasets.")


def _read_v3_10x_h5(filename: str):
    """
    Read hdf5 file from Cell Ranger v3 or later versions.
    """
    with h5py.File(str(filename), "r") as f:
        try:
            dsets = {}
            _collect_datasets(dsets, f["matrix"])

            from scipy.sparse import csr_matrix

            M, N = dsets["shape"]
            data = dsets["data"]
            if dsets["data"].dtype == np.dtype("int32"):
                data = dsets["data"].view("float32")
                data[:] = dsets["data"]
            matrix = csr_matrix(
                (data, dsets["indices"], dsets["indptr"]),
                shape=(N, M),
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                adata = anndata.AnnData(
                    matrix,
                    obs=dict(obs_names=dsets["barcodes"].astype(str)),
                    var=dict(
                        var_names=dsets["id"].astype(str),
                        gene_name=dsets["name"].astype(str),
                        feature_type=dsets["feature_type"].astype(str),
                        genome=dsets["genome"].astype(str),
                    ),
                )
            return adata
        except KeyError:
            raise Exception("File is missing one or more required datasets.")


def read_10x_h5(
    filename: t.Union[str, "pathlib.Path"],
    genome: t.Optional[str] = None,
) -> "anndata.AnnData":
    """
     Read 10x-Genomics-formatted hdf5 file.
     This method is borrowed from scanpy, please refer to scanpy library:
     https://github.com/scverse/scanpy/blob/1.6.x/scanpy/readwrite.py#LL138C10-L138C10.

    :param filename:  Path to a 10x hdf5 file.
    :param genome: Filter expression to genes within this genome. For legacy 10x h5 files, this must be provided
        if the data contains more than one genome.
    :return: Annotated data matrix, where observations/cells are named by their barcode and variables/genes by
        gene name. Stores the following information:
            :attr:`~anndata.AnnData.X`
                The data matrix is stored
            :attr:`~anndata.AnnData.obs_names`
                Cell names
            :attr:`~anndata.AnnData.var_names`
                Gene names
            :attr:`~anndata.AnnData.var`\\ `['gene_ids']`
                Gene IDs
            :attr:`~anndata.AnnData.var`\\ `['feature_types']`
                Feature types
    """
    with h5py.File(str(filename), "r") as f:
        v3 = "/matrix" in f
    if v3:
        adata = _read_v3_10x_h5(filename)
        if adata.is_view:
            adata = adata.copy()
    else:
        adata = _read_legacy_10x_h5(filename, genome=genome)
    return adata


def read_h5_or_h5ad(filename: str) -> anndata.AnnData:
    """
    Read an `h5` or `h5ad` file and return an AnnData object.

    :param filename: Path to the h5 or h5ad file.
    :return: An AnnData object.

    Usage example:

    .. code-block:: python

        import anndata

        adata = read_h5_or_h5ad("path/to/file.h5ad")
    """
    if filename.endswith(".h5"):
        return read_10x_h5(filename)
    elif filename.endswith(".h5ad"):
        return anndata.read_h5ad(filename)
    else:
        raise ValueError("File should be either a `h5` or `h5ad` file.")


def adata_to_bytes(adata: "anndata.AnnData", compression: str = "gzip") -> bytes:
    """
    Convert an :class:`anndata.AnnData` object to a byte object.

    :param adata: The AnnData object to be serialized.
    :param compression: The compression type to be used. Should be any compression type supported by h5py.
        `Default:` ``"gzip"``
    :return: Serialized bytes objects of the AnnData object.

    Usage example:

    .. code-block:: python

        import anndata

        adata = anndata.AnnData(X)
        adata_bytes = adata_to_bytes(adata)
    """

    with tempfile.NamedTemporaryFile(suffix=".h5ad") as temp_file:
        adata.write(temp_file.name, compression=compression)
        temp_file.seek(0)
        return temp_file.read()


@contextmanager
def suppress_stderr():
    original_stderr = sys.stderr
    sys.stderr = open(os.devnull, "w")
    try:
        yield
    finally:
        sys.stderr.close()
        sys.stderr = original_stderr
