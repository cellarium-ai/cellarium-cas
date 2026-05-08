"""
CLI subcommand: cellarium-cas annotate

Runs the CAS ontology-aware annotation strategy on a local .h5ad file and saves
results to an output directory.

The core logic is in :func:`annotate`, which can be imported and called directly
from Python (e.g. KFP components). The :func:`annotate_command` click command is a
thin wrapper around it.
"""

import typing as t
from pathlib import Path

import click

from cellarium.cas import constants
from cellarium.cas.cli._io import save_inferred_labels as _save_inferred_labels
from cellarium.cas.cli._io import save_metadata as _save_metadata
from cellarium.cas.cli._io import save_ontology_resource as _save_ontology_resource
from cellarium.cas.cli._io import save_ontology_response as _save_ontology_response
from cellarium.cas.client import (
    CHUNK_SIZE_ANNOTATE_DEFAULT,
    DEFAULT_PRUNE_THRESHOLD,
    DEFAULT_WEIGHTING_PREFACTOR,
    CASClient,
)

_COUNT_MATRIX_CHOICES = {v.value: v for v in constants.CountMatrixInput}


def annotate(
    input_path: str,
    output_dir: str,
    cas_api_token: str,
    cas_api_url: t.Optional[str] = None,
    cas_model_name: t.Optional[str] = None,
    chunk_size: int = CHUNK_SIZE_ANNOTATE_DEFAULT,
    count_matrix_input: str = constants.CountMatrixInput.X.value,
    feature_ids_column_name: str = "index",
    feature_names_column_name: t.Optional[str] = None,
    prune_threshold: float = DEFAULT_PRUNE_THRESHOLD,
    weighting_prefactor: float = DEFAULT_WEIGHTING_PREFACTOR,
    infer_labels: bool = False,
    min_acceptable_score: float = 0.2,
    top_k: int = 3,
    save_metadata: bool = True,
    save_ontology_resource: bool = True,
) -> t.Dict[str, t.Any]:
    """
    Annotate a single-cell dataset using the CAS ontology-aware strategy and save outputs to disk.

    This is the plain-Python entry point. The CLI command :func:`annotate_command` is a thin
    wrapper around this function and passes through all arguments unchanged.

    :param input_path: Path to the input ``.h5ad`` file.
    :param output_dir: Directory to write output files into (created if absent).
    :param cas_api_token: CAS API token.
    :param cas_api_url: CAS API URL. ``None`` uses the Cellarium production server default.
    :param cas_model_name: Model name to use for annotation. ``None`` uses the server default.
    :param chunk_size: Number of cells per annotation chunk.
    :param count_matrix_input: Where to read the count matrix from (``"X"`` or ``"raw.X"``).
    :param feature_ids_column_name: Column in ``adata.var`` containing Ensembl feature IDs, or ``"index"``.
    :param feature_names_column_name: Column in ``adata.var`` containing feature names. ``None`` skips mapping.
    :param prune_threshold: Score threshold for pruning the ontology graph output.
    :param weighting_prefactor: Weighting prefactor controlling neighbor distance weight decay.
    :param infer_labels: If ``True``, run postprocessing and save ``inferred_labels.csv``.
    :param min_acceptable_score: Minimum score for a cell type call. Used when ``infer_labels=True``.
    :param top_k: Number of top cell type calls per cell. Used when ``infer_labels=True``.
    :param save_metadata: If ``True``, save ``metadata.json`` with run provenance info.
    :param save_ontology_resource: If ``True``, save ``ontology_resource.json`` for offline benchmarking.

    :returns: Dict of paths written:
        ``output_dir``, ``ontology_response_path``, and optionally
        ``ontology_resource_path``, ``inferred_labels_path``, ``metadata_path``.
    """
    import anndata
    import pandas as pd

    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    count_matrix_enum = _COUNT_MATRIX_CHOICES[
        count_matrix_input.upper() if count_matrix_input != "raw.X" else count_matrix_input
    ]

    adata = anndata.read_h5ad(input_path)

    client_kwargs: t.Dict[str, t.Any] = {"api_token": cas_api_token}
    if cas_api_url is not None:
        client_kwargs["api_url"] = cas_api_url
    cas = CASClient(**client_kwargs)

    response = cas.annotate_matrix_cell_type_ontology_aware_strategy(
        matrix=adata,
        chunk_size=chunk_size,
        count_matrix_input=count_matrix_enum,
        feature_ids_column_name=feature_ids_column_name,
        cas_model_name=cas_model_name,
        feature_names_column_name=feature_names_column_name,
        prune_threshold=prune_threshold,
        weighting_prefactor=weighting_prefactor,
    )

    resolved_model_name = response.model_name
    result: t.Dict[str, t.Any] = {"output_dir": str(output_dir_path)}

    result["ontology_response_path"] = str(_save_ontology_response(response, output_dir_path))

    ontology_resource_name: t.Optional[str] = None
    if save_ontology_resource or infer_labels:
        ontology_resource_name = cas._resolve_ontology_resource_name(resolved_model_name)
        resource = cas.cas_api_service.get_cell_ontology_resource(ontology_resource_name)
        if save_ontology_resource:
            result["ontology_resource_path"] = str(_save_ontology_resource(resource, output_dir_path))

    if infer_labels:
        cas.insert_ontology_aware_response(response, adata)
        cas.compute_most_granular_top_k_calls_single(
            adata=adata,
            min_acceptable_score=min_acceptable_score,
            top_k=top_k,
        )

        label_cols = [f"cas_cell_type_label_{i}" for i in range(1, top_k + 1)]
        name_cols = [f"cas_cell_type_name_{i}" for i in range(1, top_k + 1)]
        score_cols = [f"cas_cell_type_score_{i}" for i in range(1, top_k + 1)]
        all_cols = [c for c in label_cols + name_cols + score_cols if c in adata.obs.columns]
        df = pd.DataFrame(adata.obs[all_cols], index=adata.obs.index)
        result["inferred_labels_path"] = str(_save_inferred_labels(df, output_dir_path))

    if save_metadata:
        metadata_dict = {
            "input_path": str(Path(input_path).resolve()),
            "model_name": resolved_model_name,
            "n_cells": len(adata),
            "ontology_resource_name": ontology_resource_name,
        }
        result["metadata_path"] = str(_save_metadata(metadata_dict, output_dir_path))

    return result


@click.command("annotate")
@click.option(
    "--input-path", required=True, type=click.Path(exists=True, dir_okay=False), help="Path to input .h5ad file."
)
@click.option(
    "--output-dir", required=True, type=click.Path(file_okay=False), help="Directory to write output files into."
)
@click.option(
    "--cas-api-token",
    default=None,
    envvar="CAS_API_TOKEN",
    help="CAS API token. Can also be set via the CAS_API_TOKEN environment variable.",
)
@click.option(
    "--cas-api-url",
    default=None,
    envvar="CAS_API_URL",
    help="CAS API URL. Defaults to the Cellarium production server. Can also be set via CAS_API_URL.",
)
@click.option(
    "--cas-model-name", default=None, help="Model name to use for annotation. Defaults to the server default."
)
@click.option(
    "--chunk-size", default=CHUNK_SIZE_ANNOTATE_DEFAULT, show_default=True, help="Number of cells per annotation chunk."
)
@click.option(
    "--count-matrix-input",
    default=constants.CountMatrixInput.X.value,
    show_default=True,
    type=click.Choice(list(_COUNT_MATRIX_CHOICES.keys()), case_sensitive=False),
    help="Where to read the count matrix from in the AnnData object.",
)
@click.option(
    "--feature-ids-column-name",
    default="index",
    show_default=True,
    help="Column in adata.var containing Ensembl feature IDs, or 'index'.",
)
@click.option(
    "--feature-names-column-name",
    default=None,
    help="Column in adata.var containing feature names (symbols). Omit to skip mapping.",
)
@click.option(
    "--prune-threshold",
    default=DEFAULT_PRUNE_THRESHOLD,
    show_default=True,
    help="Score threshold for pruning the ontology graph output.",
)
@click.option(
    "--weighting-prefactor",
    default=DEFAULT_WEIGHTING_PREFACTOR,
    show_default=True,
    help="Weighting prefactor controlling neighbor distance weight decay.",
)
@click.option(
    "--infer-labels/--no-infer-labels",
    default=False,
    show_default=True,
    help="Run postprocessing to infer top-k cell type labels and save inferred_labels.csv.",
)
@click.option(
    "--min-acceptable-score",
    default=0.2,
    show_default=True,
    help="Minimum score for a cell type call to be included. Used when --infer-labels is set.",
)
@click.option(
    "--top-k",
    default=3,
    show_default=True,
    help="Number of top cell type calls per cell. Used when --infer-labels is set.",
)
@click.option(
    "--save-metadata/--no-save-metadata",
    default=True,
    show_default=True,
    help="Save metadata.json with run provenance info.",
)
@click.option(
    "--save-ontology-resource/--no-save-ontology-resource",
    default=True,
    show_default=True,
    help="Save ontology_resource.json. Required for offline ontology-aware benchmarking.",
)
def annotate_command(
    input_path: str,
    output_dir: str,
    cas_api_token: t.Optional[str],
    cas_api_url: t.Optional[str],
    cas_model_name: t.Optional[str],
    chunk_size: int,
    count_matrix_input: str,
    feature_ids_column_name: str,
    feature_names_column_name: t.Optional[str],
    prune_threshold: float,
    weighting_prefactor: float,
    infer_labels: bool,
    min_acceptable_score: float,
    top_k: int,
    save_metadata: bool,
    save_ontology_resource: bool,
) -> None:
    """Annotate a single-cell dataset using the CAS ontology-aware strategy."""
    if cas_api_token is None:
        raise click.UsageError(
            "CAS API token is required. Pass --cas-api-token or set the CAS_API_TOKEN environment variable."
        )

    click.echo(f"Loading dataset from {input_path} ...")
    click.echo("Connecting to CAS ...")

    result = annotate(
        input_path=input_path,
        output_dir=output_dir,
        cas_api_token=cas_api_token,
        cas_api_url=cas_api_url,
        cas_model_name=cas_model_name,
        chunk_size=chunk_size,
        count_matrix_input=count_matrix_input,
        feature_ids_column_name=feature_ids_column_name,
        feature_names_column_name=feature_names_column_name,
        prune_threshold=prune_threshold,
        weighting_prefactor=weighting_prefactor,
        infer_labels=infer_labels,
        min_acceptable_score=min_acceptable_score,
        top_k=top_k,
        save_metadata=save_metadata,
        save_ontology_resource=save_ontology_resource,
    )

    click.echo(f"Saved ontology response → {result['ontology_response_path']}")
    if "ontology_resource_path" in result:
        click.echo(f"Saved ontology resource → {result['ontology_resource_path']}")
    if "inferred_labels_path" in result:
        click.echo(f"Saved inferred labels  → {result['inferred_labels_path']}")
    if "metadata_path" in result:
        click.echo(f"Saved metadata         → {result['metadata_path']}")
