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


def convert_h5ad_to_10x(input_path: str, output_dir: str) -> None:
    """
    Convert a ``.h5ad`` file to 10x Genomics sparse matrix format.

    Writes ``matrix.mtx.gz``, ``barcodes.tsv.gz``, and ``features.tsv.gz``
    into *output_dir*, suitable for loading with ``Seurat::Read10X``.

    :param input_path: Path to the input ``.h5ad`` file.
    :param output_dir: Directory to write the 10x output files into (created if absent).
    """
    os.makedirs(output_dir, exist_ok=True)

    adata = anndata.read_h5ad(input_path)

    scipy.io.mmwrite(f"{output_dir}/matrix.mtx", adata.X.T.tocsc())
    with open(f"{output_dir}/matrix.mtx", "rb") as f_in, gzip.open(f"{output_dir}/matrix.mtx.gz", "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(f"{output_dir}/matrix.mtx")

    adata.obs_names.to_series().to_csv(f"{output_dir}/barcodes.tsv.gz", index=False, header=False, compression="gzip")

    features = pd.DataFrame({"gene_id": adata.var_names, "gene_name": adata.var_names})
    features.to_csv(f"{output_dir}/features.tsv.gz", index=False, header=False, compression="gzip")

    print("Conversion done.")


if __name__ == "__main__":
    convert_h5ad_to_10x(input_path=sys.argv[1], output_dir=sys.argv[2])
