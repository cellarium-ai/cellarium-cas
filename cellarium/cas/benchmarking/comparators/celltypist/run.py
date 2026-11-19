#!/usr/bin/env python3
"""Run raw per-cell CellTypist inference and write native ranked predictions."""

from __future__ import annotations

import hashlib
import json
from importlib.metadata import version
from pathlib import Path
from typing import Any, Sequence

import anndata
import click
import numpy as np
import pandas as pd
from scipy import sparse

PREDICTIONS_FILENAME = "celltypist_predictions.csv"
STATS_FILENAME = "celltypist_stats.json"
LABELS_FILENAME = "predicted_labels.csv"
MODE = "best match"
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
    """Build log1p-normalized (10,000 counts per cell), symbol-consolidated AnnData for CellTypist.

    Raw counts of repeated symbols are summed, rows are scaled across all input genes, and log1p
    is applied. Zero-total cells stay all-zero.

    :param adata: Read raw counts, barcodes, and gene symbols from this AnnData object.

    :raises ValueError: If cells, genes, barcodes, symbols, or counts are invalid.
    :returns: Return AnnData with barcode index and symbol ``var_names``.
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
        var=pd.DataFrame(index=pd.Index(unique_symbols.astype(str))),
    )


def _rank_native(
    *,
    predicted_labels: pd.DataFrame,
    decision: pd.DataFrame,
    probability: pd.DataFrame,
    barcodes: pd.Index,
    top_k: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Align upstream results to input order and rank native classes by decision score.

    Ties break by model class-column order. Scores written are the original sigmoid probabilities;
    ranks beyond the class count are empty.

    :param predicted_labels: Read upstream ``predicted_labels`` table (column ``predicted_labels``).
    :param decision: Read upstream cells x classes decision matrix.
    :param probability: Read upstream cells x classes probability matrix.
    :param barcodes: Align rows to these input barcodes.
    :param top_k: Write this many ranks.

    :raises ValueError: If results are misaligned, non-finite, out of range, or disagree with upstream labels.
    :returns: Return ``(aligned predicted_labels, ranked predictions)``.
    """
    aligned: dict[str, pd.DataFrame] = {}
    for name, frame in (
        ("predicted_labels", predicted_labels),
        ("decision_matrix", decision),
        ("probability_matrix", probability),
    ):
        index = pd.Index(frame.index.astype(str))
        if not index.is_unique or len(index) != len(barcodes) or set(index) != set(barcodes):
            raise ValueError(f"CellTypist {name} barcodes do not match the input barcodes one-to-one.")
        frame = frame.copy()
        frame.index = index
        aligned[name] = frame.reindex(barcodes)
    labels_df, decision, probability = (
        aligned["predicted_labels"],
        aligned["decision_matrix"],
        aligned["probability_matrix"],
    )

    if list(decision.columns) != list(probability.columns):
        raise ValueError("CellTypist decision and probability matrices have different class columns.")
    scores = decision.to_numpy(dtype=np.float64)
    probs = probability.to_numpy(dtype=np.float64)
    if not np.isfinite(scores).all() or not np.isfinite(probs).all():
        raise ValueError("CellTypist returned non-finite decision scores or probabilities.")
    if (probs < 0).any() or (probs > 1).any():
        raise ValueError("CellTypist probabilities fall outside [0, 1].")

    classes = np.asarray(decision.columns.astype(str), dtype=object)
    order = np.argsort(-scores, axis=1, kind="stable")[:, :top_k]
    if (
        "predicted_labels" not in labels_df.columns
        or not (labels_df["predicted_labels"].astype(str).to_numpy() == classes[order[:, 0]]).all()
    ):
        raise ValueError("Decision-score top-1 disagrees with CellTypist predicted_labels.")

    ranked = pd.DataFrame(index=barcodes)
    rows = np.arange(len(barcodes))[:, None]
    top_labels = classes[order]
    top_probs = probs[rows, order]
    for rank in range(top_k):
        if rank < order.shape[1]:
            ranked[f"predicted_label_{rank + 1}"] = top_labels[:, rank]
            ranked[f"predicted_score_{rank + 1}"] = top_probs[:, rank]
        else:
            ranked[f"predicted_label_{rank + 1}"] = np.nan
            ranked[f"predicted_score_{rank + 1}"] = np.nan
    return labels_df, ranked


def _sha256(path: Path) -> str:
    """Return the hex SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with open(path, "rb") as file:
        for chunk in iter(lambda: file.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def _write_outputs(*, output_dir: Path, labels: pd.DataFrame, ranked: pd.DataFrame, stats: dict[str, Any]) -> None:
    """Write raw upstream labels, ranked predictions, and statistics JSON."""
    output_dir.mkdir(parents=True, exist_ok=True)
    labels.to_csv(output_dir / LABELS_FILENAME)
    ranked.to_csv(output_dir / PREDICTIONS_FILENAME)
    with open(output_dir / STATS_FILENAME, "w") as file:
        json.dump(stats, file, indent=2, default=str)


def run_celltypist(
    *,
    model_path: str,
    input_h5ad: tuple[str, ...],
    output_dir_template: str,
    top_k: int = 5,
) -> list[dict[str, object]]:
    """Run raw per-cell CellTypist (``best match``, no majority voting) over one or more H5AD inputs.

    Only trusted model files may be used: loading a model unpickles it.

    :param model_path: Load this existing CellTypist ``.pkl`` model file.
    :param input_h5ad: Annotate these raw-count H5AD files.
    :param output_dir_template: Write each result under this fresh directory template containing ``{stem}``.
    :param top_k: Write this many ranked native classes per cell.

    :raises ValueError: If options, inputs, counts, gene overlap, or native results are invalid.
        Earlier completed inputs remain written.
    :returns: Return per-input statistics.
    """
    if top_k <= 0:
        raise ValueError("top_k must be positive.")
    model_file = Path(model_path).resolve()
    if not model_file.is_file():
        raise ValueError(f"Model file does not exist: {model_path}")
    output_dirs = _resolve_output_dirs(input_h5ad=input_h5ad, output_dir_template=output_dir_template)

    import celltypist  # deferred: not a base CAS dependency

    model = celltypist.models.Model.load(model=str(model_file))
    model_sha256 = _sha256(model_file)
    model_genes = set(map(str, model.features))
    api_version = version("celltypist")
    model_name = f"celltypist_{model_file.stem}_{model_sha256}_raw"

    stats_list: list[dict[str, object]] = []
    for input_path, output_dir in zip(input_h5ad, output_dirs):
        adata = anndata.read_h5ad(filename=input_path)
        n_input_genes = int(adata.n_vars)
        prepared = _prepare_dataset(adata=adata)
        del adata

        matched = int(prepared.var_names.isin(model_genes).sum())
        if matched == 0:
            raise ValueError(
                f"No gene symbols in '{input_path}' overlap the CellTypist model features. "
                "Check that gene symbols are used."
            )

        result = celltypist.annotate(filename=prepared, model=model, mode=MODE, majority_voting=False)
        labels, ranked = _rank_native(
            predicted_labels=result.predicted_labels,
            decision=result.decision_matrix,
            probability=result.probability_matrix,
            barcodes=prepared.obs_names,
            top_k=top_k,
        )
        stats: dict[str, object] = {
            "input_h5ad": str(Path(input_path).resolve()),
            "n_cells": int(prepared.n_obs),
            "n_input_genes": n_input_genes,
            "n_model_genes": len(model_genes),
            "n_matched_genes": matched,
            "api_version": api_version,
            "model_path": str(model_file),
            "model_sha256": model_sha256,
            "model_description": getattr(model, "description", None),
            "model_name": model_name,
            "mode": MODE,
            "majority_voting": False,
            "top_k": top_k,
        }
        _write_outputs(output_dir=output_dir, labels=labels, ranked=ranked, stats=stats)
        stats_list.append({"stem": Path(input_path).stem, **stats})
    return stats_list


@click.command("run-celltypist")
@click.option("--model-path", required=True, type=click.Path(exists=True, dir_okay=False, resolve_path=True))
@click.option(
    "--input-h5ad", required=True, multiple=True, type=click.Path(exists=True, dir_okay=False, resolve_path=True)
)
@click.option("--output-dir-template", required=True)
@click.option("--top-k", default=5, show_default=True, type=click.IntRange(min=1))
def run_celltypist_command(model_path: str, input_h5ad: tuple[str, ...], output_dir_template: str, top_k: int) -> None:
    """Annotate one or more raw-count H5AD files with a CellTypist model (best match, no majority voting)."""
    try:
        stats_list = run_celltypist(
            model_path=model_path, input_h5ad=input_h5ad, output_dir_template=output_dir_template, top_k=top_k
        )
    except ValueError as error:
        raise click.ClickException(str(error)) from error
    for stats in stats_list:
        click.echo(
            f"[celltypist] {stats['stem']}: cells={stats['n_cells']} "
            f"matched_genes={stats['n_matched_genes']}/{stats['n_model_genes']}"
        )


if __name__ == "__main__":
    run_celltypist_command()
