#!/usr/bin/env python3
"""Run SCimilarity inference and write ranked predictions."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Sequence

import anndata
import click
import numpy as np
import pandas as pd
from scipy import sparse


def _load_scimilarity_utilities() -> tuple[Any, Any, Any]:
    """Load SCimilarity preprocessing utilities at runtime."""
    from scimilarity.utils import align_dataset, consolidate_duplicate_symbols, lognorm_counts

    return align_dataset, consolidate_duplicate_symbols, lognorm_counts


def _gene_symbols(*, adata: anndata.AnnData) -> pd.Index:
    """Select gene symbols using the documented column precedence."""
    for column in ("feature_name", "gene_symbol"):
        if column in adata.var.columns:
            symbols = adata.var[column]
            if symbols.isna().any() or (symbols.astype(str).str.strip() == "").any():
                raise ValueError(f"Gene symbol column '{column}' contains missing or empty values.")
            return pd.Index(symbols.astype(str))
    return pd.Index(adata.var_names.astype(str))


def _raw_counts(*, adata: anndata.AnnData) -> Any:
    """Return raw counts, preferring the counts layer."""
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if counts is None:
        raise ValueError("Raw counts are required in layers['counts'] or X.")
    if sparse.issparse(counts):
        if counts.data.size and np.min(counts.data) < 0:
            raise ValueError("Raw counts must not contain negative values.")
    elif np.any(np.asarray(counts) < 0):
        raise ValueError("Raw counts must not contain negative values.")
    return counts


def _prepare_dataset(
    *,
    adata: anndata.AnnData,
    gene_order: Sequence[str],
    gene_overlap_threshold: int,
) -> tuple[anndata.AnnData, int, int, int]:
    """Prepare raw-count data for SCimilarity embedding.

    :param adata: Read raw-count observations from this AnnData object.
    :param gene_order: Align genes to this ordered model gene list.
    :param gene_overlap_threshold: Require at least this many matched genes.

    :raises ValueError: If input has no cells, raw counts, or sufficient model genes.
    :returns: Return normalized data and input, model, and matched gene counts.
    """
    if adata.n_obs == 0:
        raise ValueError("Input H5AD contains no cells.")

    _raw_counts(adata=adata)
    input_gene_count = int(adata.n_vars)
    data = adata
    if "counts" not in data.layers:
        data.layers["counts"] = data.X
    data.X = data.layers["counts"]
    data.var_names = _gene_symbols(adata=adata)

    align_dataset, consolidate_duplicate_symbols, lognorm_counts = _load_scimilarity_utilities()
    data = consolidate_duplicate_symbols(adata=data)
    matched_gene_count = int(data.var_names.isin(gene_order).sum())
    if matched_gene_count < gene_overlap_threshold:
        raise ValueError(
            f"Dataset incompatible. Gene overlap of {matched_gene_count} is less than "
            f"--gene-overlap-threshold {gene_overlap_threshold}. Check that gene symbols are used."
        )
    data = align_dataset(
        data=data,
        target_gene_order=gene_order,
        gene_overlap_threshold=gene_overlap_threshold,
    )
    data = lognorm_counts(data=data)
    return data, input_gene_count, len(gene_order), matched_gene_count


def _rank_weighted_hits(*, weighted_hits: str, top_k: int) -> list[tuple[str, float]]:
    """Rank normalized inverse-distance votes with stable tie ordering."""
    try:
        votes = json.loads(weighted_hits)
    except (TypeError, json.JSONDecodeError) as error:
        raise ValueError("Malformed SCimilarity hits_weighted data.") from error
    if not isinstance(votes, dict) or not votes:
        raise ValueError("Malformed SCimilarity hits_weighted data.")

    parsed_votes: list[tuple[str, float]] = []
    for label, vote in votes.items():
        if not isinstance(label, str) or not isinstance(vote, (int, float)) or not math.isfinite(vote) or vote < 0:
            raise ValueError("Malformed SCimilarity hits_weighted data.")
        parsed_votes.append((label, float(vote)))
    total_vote = sum(vote for _, vote in parsed_votes)
    if total_vote <= 0:
        raise ValueError("Malformed SCimilarity hits_weighted data.")
    return [
        (label, vote / total_vote)
        for label, vote in sorted(parsed_votes, key=lambda item: item[1], reverse=True)[:top_k]
    ]


def _predictions_dataframe(
    *,
    obs_names: Sequence[str],
    predictions: pd.Series,
    prediction_stats: pd.DataFrame,
    top_k: int,
) -> pd.DataFrame:
    """Build top-k prediction columns from SCimilarity weighted neighbor hits."""
    if len(predictions) != len(obs_names) or len(prediction_stats) != len(obs_names):
        raise ValueError("SCimilarity predictions do not match input cell count.")
    if "hits_weighted" not in prediction_stats:
        raise ValueError("SCimilarity prediction statistics lack hits_weighted.")

    columns = [
        column for rank in range(1, top_k + 1) for column in (f"predicted_label_{rank}", f"predicted_score_{rank}")
    ]
    rows: list[list[str | float | None]] = []
    for position, weighted_hits in enumerate(prediction_stats["hits_weighted"]):
        ranked = _rank_weighted_hits(weighted_hits=weighted_hits, top_k=top_k)
        if str(predictions.iloc[position]) != ranked[0][0]:
            raise ValueError("SCimilarity weighted prediction disagrees with hits_weighted ranking.")
        row: list[str | float | None] = []
        for rank in range(top_k):
            if rank < len(ranked):
                row.extend(ranked[rank])
            else:
                row.extend((None, None))
        rows.append(row)
    return pd.DataFrame(rows, index=pd.Index(obs_names, name=None), columns=columns)


def _validate_options(
    *,
    input_h5ad: Sequence[str],
    output_dir_template: str,
    top_k: int,
    neighbors: int,
    ef: int,
    gene_overlap_threshold: int,
    buffer_size: int,
) -> None:
    """Reject invalid runner options before loading the model."""
    if not input_h5ad:
        raise ValueError("At least one input H5AD is required.")
    if "{stem}" not in output_dir_template:
        raise ValueError("--output-dir-template must contain '{stem}'.")
    for option_name, value in {
        "top_k": top_k,
        "neighbors": neighbors,
        "ef": ef,
        "gene_overlap_threshold": gene_overlap_threshold,
        "buffer_size": buffer_size,
    }.items():
        if value <= 0:
            raise ValueError(f"{option_name} must be positive.")


def run_scimilarity(
    *,
    model_path: str,
    input_h5ad: Sequence[str],
    output_dir_template: str,
    top_k: int = 5,
    neighbors: int = 50,
    ef: int = 100,
    gene_overlap_threshold: int = 5000,
    buffer_size: int = 10000,
) -> list[dict[str, Any]]:
    """Run unconstrained SCimilarity inference over one or more H5AD inputs.

    :param model_path: Load the SCimilarity model from this directory.
    :param input_h5ad: Process these H5AD files with one loaded model.
    :param output_dir_template: Write each result under this template containing ``{stem}``.
    :param top_k: Write this many weighted-neighbor labels per cell.
    :param neighbors: Query this many nearest neighbors.
    :param ef: Use this HNSW dynamic-list size.
    :param gene_overlap_threshold: Require this many matched model genes.
    :param buffer_size: Embed this many cells per batch.

    :raises ValueError: If options, input counts, gene overlap, or predictions are invalid.
    :returns: Return per-input inference statistics.
    """
    _validate_options(
        input_h5ad=input_h5ad,
        output_dir_template=output_dir_template,
        top_k=top_k,
        neighbors=neighbors,
        ef=ef,
        gene_overlap_threshold=gene_overlap_threshold,
        buffer_size=buffer_size,
    )

    from scimilarity import CellAnnotation

    ca = CellAnnotation(model_path=model_path, use_gpu=False)
    if top_k > len(ca.classes):
        raise ValueError(f"top_k {top_k} exceeds the model's {len(ca.classes)} classes.")

    stats_list: list[dict[str, Any]] = []
    for input_path in input_h5ad:
        adata = anndata.read_h5ad(filename=input_path)
        obs_names = adata.obs_names.astype(str).tolist()
        n_cells = int(adata.n_obs)
        data, input_gene_count, model_gene_count, matched_gene_count = _prepare_dataset(
            adata=adata,
            gene_order=ca.gene_order,
            gene_overlap_threshold=gene_overlap_threshold,
        )
        embeddings = ca.get_embeddings(X=data.X, buffer_size=buffer_size)
        predictions, _, _, prediction_stats = ca.get_predictions_knn(
            embeddings=embeddings,
            k=neighbors,
            ef=ef,
            weighting=True,
            disable_progress=True,
        )
        predictions_df = _predictions_dataframe(
            obs_names=obs_names,
            predictions=predictions,
            prediction_stats=prediction_stats,
            top_k=top_k,
        )

        output_dir = Path(output_dir_template.format(stem=Path(input_path).stem))
        output_dir.mkdir(parents=True, exist_ok=True)
        predictions_df.to_csv(output_dir / "scimilarity_predictions.csv")
        stats = {
            "input_h5ad": str(Path(input_path).resolve()),
            "n_cells": n_cells,
            "n_input_genes": input_gene_count,
            "n_model_genes": model_gene_count,
            "n_matched_genes": matched_gene_count,
            "overlap_ratio": float(matched_gene_count / model_gene_count),
            "api_version": "0.4.1",
            "model_version": "model_v1.1",
            "neighbors": neighbors,
            "ef": ef,
            "weighting": "inverse_distance",
        }
        with open(output_dir / "scimilarity_stats.json", "w") as file:
            json.dump(stats, file, indent=2)
        stats_list.append({"stem": Path(input_path).stem, **stats})
    return stats_list


@click.command("run-scimilarity")
@click.option("--model-path", required=True, type=click.Path(exists=True, file_okay=False, resolve_path=True))
@click.option(
    "--input-h5ad", required=True, multiple=True, type=click.Path(exists=True, dir_okay=False, resolve_path=True)
)
@click.option("--output-dir-template", required=True)
@click.option("--top-k", default=5, show_default=True, type=click.IntRange(min=1))
@click.option("--neighbors", default=50, show_default=True, type=click.IntRange(min=1))
@click.option("--ef", default=100, show_default=True, type=click.IntRange(min=1))
@click.option("--gene-overlap-threshold", default=5000, show_default=True, type=click.IntRange(min=1))
@click.option("--buffer-size", default=10000, show_default=True, type=click.IntRange(min=1))
def run_scimilarity_command(
    model_path: str,
    input_h5ad: Sequence[str],
    output_dir_template: str,
    top_k: int,
    neighbors: int,
    ef: int,
    gene_overlap_threshold: int,
    buffer_size: int,
) -> None:
    """Run SCimilarity inference over one or more H5AD inputs."""
    stats_list = run_scimilarity(
        model_path=model_path,
        input_h5ad=input_h5ad,
        output_dir_template=output_dir_template,
        top_k=top_k,
        neighbors=neighbors,
        ef=ef,
        gene_overlap_threshold=gene_overlap_threshold,
        buffer_size=buffer_size,
    )
    for stats in stats_list:
        click.echo(
            f"[scimilarity] {stats['stem']}: cells={stats['n_cells']} "
            f"matched_genes={stats['n_matched_genes']}/{stats['n_model_genes']}"
        )


if __name__ == "__main__":
    run_scimilarity_command()
