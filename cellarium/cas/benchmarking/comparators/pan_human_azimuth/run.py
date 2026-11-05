#!/usr/bin/env python3
"""Run Pan-Human Azimuth inference and write official native annotations."""

from __future__ import annotations

import json
from importlib.metadata import version
from importlib.resources import read_text
from pathlib import Path
from typing import Any, Sequence

import anndata
import click
import numpy as np
import pandas as pd
from scipy import sparse

MODEL_VERSION = "v1"
PREDICTIONS_FILENAME = "pan_human_azimuth_predictions.csv"
STATS_FILENAME = "pan_human_azimuth_stats.json"
LEVEL_COLUMNS = ("azimuth_fine", "azimuth_medium", "azimuth_broad")
NATIVE_COLUMNS = (
    *LEVEL_COLUMNS,
    "full_hierarchical_labels",
    "final_level_labels",
    "final_level_confidence",
    "full_consistent_hierarchy",
)
UNASSIGNED = "Unassigned"
MAX_SAMPLE_BARCODES = 5
TARGET_SUM = 10_000


def _barcodes(*, adata: anndata.AnnData) -> pd.Index:
    """Return validated observation identifiers.

    :param adata: Read observation identifiers from this AnnData object.

    :raises ValueError: If there are no cells or identifiers are blank or repeated.
    :returns: Return string barcodes in input order.
    """
    if adata.n_obs == 0:
        raise ValueError("Input H5AD contains no cells.")
    if adata.n_vars == 0:
        raise ValueError("Input H5AD contains no genes.")
    barcodes = pd.Index(adata.obs_names.astype(str))
    if (barcodes.str.strip() == "").any():
        raise ValueError("Input H5AD contains blank cell barcodes.")
    if not barcodes.is_unique:
        raise ValueError("Input H5AD contains duplicate cell barcodes.")
    return barcodes


def _gene_symbols(*, adata: anndata.AnnData) -> pd.Index:
    """Select gene symbols using ``feature_name``, ``gene_symbol``, then ``var_names``."""
    for column in ("feature_name", "gene_symbol"):
        if column in adata.var.columns:
            symbols = adata.var[column]
            if symbols.isna().any() or (symbols.astype(str).str.strip() == "").any():
                raise ValueError(f"Gene symbol column '{column}' contains missing or empty values.")
            return pd.Index(symbols.astype(str))
    symbols = pd.Index(adata.var_names.astype(str))
    if (symbols.str.strip() == "").any():
        raise ValueError("var_names contains empty gene symbols.")
    return symbols


def _raw_counts(*, adata: anndata.AnnData) -> sparse.csr_matrix:
    """Return validated raw counts from ``layers['counts']``, else ``X``, as CSR.

    :raises ValueError: If counts are absent, non-finite, negative, or fractional.
    """
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if counts is None:
        raise ValueError("Raw counts are required in layers['counts'] or X.")
    counts = sparse.csr_matrix(counts)
    values = counts.data
    if not np.isfinite(values).all():
        raise ValueError("Raw counts must be finite.")
    if (values < 0).any():
        raise ValueError("Raw counts must not contain negative values.")
    if (values != np.floor(values)).any():
        raise ValueError("Raw counts must be integer-valued; got fractional (normalized?) values.")
    return counts


def _prepare_dataset(*, adata: anndata.AnnData) -> anndata.AnnData:
    """Build log-normalized, symbol-consolidated AnnData for Pan-Human Azimuth.

    Raw counts of repeated symbols are summed, rows are scaled to 10,000 counts across all input
    genes, and log1p is applied. Zero-total cells stay all-zero.

    :param adata: Read raw counts, barcodes, and gene symbols from this AnnData object.

    :raises ValueError: If cells, genes, barcodes, symbols, or counts are invalid.
    :returns: Return AnnData with barcode index and ``var['feature_name']`` symbols.
    """
    barcodes = _barcodes(adata=adata)
    counts = _raw_counts(adata=adata)
    symbols = _gene_symbols(adata=adata)

    codes, unique_symbols = pd.factorize(symbols)  # first-occurrence order
    if len(unique_symbols) < len(symbols):
        aggregation = sparse.csr_matrix(
            (np.ones(len(codes)), (np.arange(len(codes)), codes)), shape=(len(codes), len(unique_symbols))
        )
        counts = (counts @ aggregation).tocsr()

    totals = np.asarray(counts.sum(axis=1), dtype=np.float64).ravel()
    scale = np.divide(TARGET_SUM, totals, out=np.zeros_like(totals), where=totals > 0)
    normalized = (sparse.diags(scale) @ counts).tocsr()
    normalized.data = np.log1p(normalized.data)

    return anndata.AnnData(
        X=normalized.astype(np.float32),
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame({"feature_name": unique_symbols.astype(str)}, index=pd.Index(unique_symbols.astype(str))),
    )


def _model_genes() -> set[str]:
    """Read the packaged v1 feature panel (importing ``panhumanpy`` loads TensorFlow)."""
    text = read_text(f"panhumanpy._tools.{MODEL_VERSION}.inference_feature_panel", "inference_feature_panel.txt")
    return {line.strip() for line in text.splitlines() if line.strip()}


def _api_version() -> str:
    """Return installed ``panhumanpy`` package version."""
    return version("panhumanpy")


def _annotate(*, prepared: anndata.AnnData, batch_size: int) -> pd.DataFrame:
    """Run official ``AzimuthNN`` and return its native per-cell metadata.

    This is the inference boundary.

    # ponytail: high-level API reloads model objects per input; use AzimuthNN_base to retain
    # them across files if repeated-file startup becomes a measured problem.
    """
    from panhumanpy import AzimuthNN

    azimuth = AzimuthNN(
        query_arg=prepared,
        feature_names_col="feature_name",
        model_version=MODEL_VERSION,
        eval_batch_size=batch_size,
        normalization_override=True,
        output_mode="minimal",
        refine=True,
    )
    return azimuth.cells_meta


def _validate_native(*, cells_meta: pd.DataFrame, barcodes: pd.Index, input_h5ad: str) -> pd.DataFrame:
    """Align native metadata to input barcodes and reject invalid or ``Unassigned`` annotations.

    :param cells_meta: Native metadata returned by Pan-Human Azimuth.
    :param barcodes: Input barcodes in original order.
    :param input_h5ad: Name the input in error messages.

    :raises ValueError: If barcodes or native columns are invalid, or any level is ``Unassigned``.
    :returns: Return the seven native columns indexed by ``barcodes``.
    """
    index = pd.Index(cells_meta.index.astype(str))
    if not index.is_unique:
        raise ValueError("Native annotations contain duplicate barcodes.")
    if set(index) != set(barcodes):
        missing = len(set(barcodes) - set(index))
        extra = len(set(index) - set(barcodes))
        raise ValueError(f"Native annotation barcodes differ from input: {missing} missing, {extra} unexpected.")
    absent = [column for column in NATIVE_COLUMNS if column not in cells_meta.columns]
    if absent:
        raise ValueError(f"Native annotations are missing columns: {absent}.")

    native = cells_meta.loc[:, list(NATIVE_COLUMNS)].copy()
    native.index = index
    native = native.reindex(barcodes)
    for column in LEVEL_COLUMNS:
        labels = native[column]
        if labels.isna().any() or (labels.astype(str).str.strip() == "").any():
            raise ValueError(f"Native column '{column}' contains missing or empty labels.")

    unassigned = np.zeros(len(native), dtype=bool)
    for column in LEVEL_COLUMNS:
        unassigned |= (native[column].astype(str) == UNASSIGNED).to_numpy()
    if unassigned.any():
        affected = barcodes[unassigned]
        sample = list(affected[:MAX_SAMPLE_BARCODES])
        raise ValueError(
            f"Pan-Human Azimuth produced Unassigned annotations for {len(affected)} cells in "
            f"'{input_h5ad}' (first {len(sample)} barcodes): {sample}"
        )
    return native


def _resolve_output_dirs(*, input_h5ad: Sequence[str], output_dir_template: str) -> list[Path]:
    """Validate options and return resolved, fresh output directories, one per input.

    :raises ValueError: If options are invalid, inputs are missing, or outputs collide or exist nonempty.
    """
    if not input_h5ad:
        raise ValueError("At least one input H5AD is required.")
    if "{stem}" not in output_dir_template:
        raise ValueError("--output-dir-template must contain '{stem}'.")
    output_dirs: list[Path] = []
    for input_path in input_h5ad:
        if not Path(input_path).is_file():
            raise ValueError(f"Input H5AD does not exist: {input_path}")
        output_dirs.append(Path(output_dir_template.format(stem=Path(input_path).stem)).resolve())
    if len(set(output_dirs)) != len(output_dirs):
        raise ValueError("Multiple inputs resolve to the same output directory.")
    for output_dir in output_dirs:
        if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
            raise ValueError(f"Output directory exists and is not empty: {output_dir}")
    return output_dirs


def _write_outputs(*, output_dir: Path, native: pd.DataFrame, stats: dict[str, Any]) -> None:
    """Write native predictions CSV and statistics JSON."""
    output_dir.mkdir(parents=True, exist_ok=True)
    native.to_csv(output_dir / PREDICTIONS_FILENAME)
    with open(output_dir / STATS_FILENAME, "w") as file:
        json.dump(stats, file, indent=2)


def run_pan_human_azimuth(
    *,
    input_h5ad: tuple[str, ...],
    output_dir_template: str,
    batch_size: int = 8192,
) -> list[dict[str, object]]:
    """Run official Pan-Human Azimuth v1 over one or more H5AD inputs.

    :param input_h5ad: Annotate these raw-count H5AD files.
    :param output_dir_template: Write each result under this fresh directory template containing ``{stem}``.
    :param batch_size: Evaluate this many cells per inference minibatch.

    :raises ValueError: If options, inputs, counts, gene overlap, or native annotations are invalid,
        or any native annotation is ``Unassigned``. Earlier completed inputs remain written.
    :returns: Return per-input statistics.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be positive.")
    output_dirs = _resolve_output_dirs(input_h5ad=input_h5ad, output_dir_template=output_dir_template)

    stats_list: list[dict[str, object]] = []
    for input_path, output_dir in zip(input_h5ad, output_dirs):
        adata = anndata.read_h5ad(filename=input_path)
        n_input_genes = int(adata.n_vars)
        prepared = _prepare_dataset(adata=adata)
        del adata

        model_genes = _model_genes()
        matched = int(prepared.var["feature_name"].isin(model_genes).sum())
        if matched == 0:
            raise ValueError(
                f"No gene symbols in '{input_path}' overlap the Pan-Human Azimuth feature panel. "
                "Check that gene symbols are used."
            )

        cells_meta = _annotate(prepared=prepared, batch_size=batch_size)
        native = _validate_native(
            cells_meta=cells_meta, barcodes=prepared.obs_names, input_h5ad=str(Path(input_path).resolve())
        )
        stats: dict[str, object] = {
            "input_h5ad": str(Path(input_path).resolve()),
            "n_cells": int(prepared.n_obs),
            "n_input_genes": n_input_genes,
            "n_model_genes": len(model_genes),
            "n_matched_genes": matched,
            "overlap_ratio": matched / len(model_genes),
            "api_version": _api_version(),
            "model_version": MODEL_VERSION,
        }
        _write_outputs(output_dir=output_dir, native=native, stats=stats)
        stats_list.append({"stem": Path(input_path).stem, **stats})
    return stats_list


@click.command("run-pan-human-azimuth")
@click.option(
    "--input-h5ad", required=True, multiple=True, type=click.Path(exists=True, dir_okay=False, resolve_path=True)
)
@click.option("--output-dir-template", required=True)
@click.option("--batch-size", default=8192, show_default=True, type=click.IntRange(min=1))
def run_pan_human_azimuth_command(input_h5ad: tuple[str, ...], output_dir_template: str, batch_size: int) -> None:
    """Annotate one or more raw-count H5AD files with Pan-Human Azimuth v1."""
    try:
        stats_list = run_pan_human_azimuth(
            input_h5ad=input_h5ad, output_dir_template=output_dir_template, batch_size=batch_size
        )
    except ValueError as error:
        raise click.ClickException(str(error)) from error
    for stats in stats_list:
        click.echo(
            f"[pan_human_azimuth] {stats['stem']}: cells={stats['n_cells']} "
            f"matched_genes={stats['n_matched_genes']}/{stats['n_model_genes']}"
        )


if __name__ == "__main__":
    run_pan_human_azimuth_command()
