#!/usr/bin/env python3
"""
Convert a .h5ad file to 10x Genomics sparse matrix format (matrix.mtx.gz,
barcodes.tsv.gz, features.tsv.gz) for consumption by Seurat / Azimuth.

Usage:
    python cellarium/cas/benchmarking/azimuth/helpers/convert_h5ad.py <input_h5ad_path> <output_dir>
"""

import gzip
import os
import shutil
import sys
from pathlib import Path

import anndata
import pandas as pd
import scipy.io
from scipy import sparse

from cellarium.cas._path_utils import resolve_within as _resolve_within


def _as_clean_str_series(values) -> pd.Series:
    s = pd.Series(values).astype("string").fillna("")
    s = s.replace({"nan": "", "None": "", "<NA>": ""})
    return s.astype(str)


def convert_h5ad_to_10x(input_path: str, output_dir: str) -> None:
    input_path = Path(input_path).resolve()
    output_dir_path = Path(output_dir).resolve()
    output_dir_path.mkdir(parents=True, exist_ok=True)

    adata = anndata.read_h5ad(input_path)

    X = adata.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)

    mtx_file = _resolve_within(output_dir_path, "matrix.mtx")
    mtx_gz_file = _resolve_within(output_dir_path, "matrix.mtx.gz")
    scipy.io.mmwrite(str(mtx_file), X.T.tocsc())

    with open(mtx_file, "rb") as f_in, gzip.open(mtx_gz_file, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    os.remove(mtx_file)

    adata.obs_names.to_series().astype(str).to_csv(
        _resolve_within(output_dir_path, "barcodes.tsv.gz"),
        index=False,
        header=False,
        compression="gzip",
    )

    if "gene_id" in adata.var.columns:
        gene_ids = _as_clean_str_series(adata.var["gene_id"])
    elif "feature_id" in adata.var.columns:
        gene_ids = _as_clean_str_series(adata.var["feature_id"])
    elif "gene_ids" in adata.var.columns:
        gene_ids = _as_clean_str_series(adata.var["gene_ids"])
    else:
        gene_ids = _as_clean_str_series(adata.var_names)

    if "feature_name" in adata.var.columns:
        gene_names = _as_clean_str_series(adata.var["feature_name"])
    elif "gene_symbol" in adata.var.columns:
        gene_names = _as_clean_str_series(adata.var["gene_symbol"])
    elif "gene_symbols" in adata.var.columns:
        gene_names = _as_clean_str_series(adata.var["gene_symbols"])
    else:
        gene_names = _as_clean_str_series(adata.var_names)

    # fallback if some IDs are empty
    gene_ids = gene_ids.mask(gene_ids == "", gene_names)

    features = pd.DataFrame(
        {
            "gene_id": gene_ids.to_numpy(),
            "gene_name": gene_names.to_numpy(),
        }
    )

    features.to_csv(
        _resolve_within(output_dir_path, "features.tsv.gz"),
        sep="\t",
        index=False,
        header=False,
        compression="gzip",
    )

    print("Conversion done.")


if __name__ == "__main__":
    convert_h5ad_to_10x(input_path=sys.argv[1], output_dir=sys.argv[2])
