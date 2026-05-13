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

import anndata
import pandas as pd
import scipy.io
from scipy import sparse


def _as_clean_str_series(values) -> pd.Series:
    s = pd.Series(values).astype("string").fillna("")
    s = s.replace({"nan": "", "None": "", "<NA>": ""})
    return s.astype(str)


def convert_h5ad_to_10x(input_path: str, output_dir: str) -> None:
    os.makedirs(output_dir, exist_ok=True)

    adata = anndata.read_h5ad(input_path)

    X = adata.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)

    scipy.io.mmwrite(f"{output_dir}/matrix.mtx", X.T.tocsc())

    with open(f"{output_dir}/matrix.mtx", "rb") as f_in, gzip.open(f"{output_dir}/matrix.mtx.gz", "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    os.remove(f"{output_dir}/matrix.mtx")

    adata.obs_names.to_series().astype(str).to_csv(
        f"{output_dir}/barcodes.tsv.gz",
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
        f"{output_dir}/features.tsv.gz",
        sep="\t",
        index=False,
        header=False,
        compression="gzip",
    )

    print("Conversion done.")


if __name__ == "__main__":
    convert_h5ad_to_10x(input_path=sys.argv[1], output_dir=sys.argv[2])
