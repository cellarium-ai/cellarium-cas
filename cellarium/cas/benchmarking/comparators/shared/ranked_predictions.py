"""Map ranked comparator predictions to CAS label files."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import anndata
import pandas as pd

from cellarium.cas.benchmarking.comparators.shared.crosswalk import load_crosswalk_maps
from cellarium.cas.postprocessing.cell_ontology.cell_ontology_cache import CellOntologyCache


def _load_and_align_predictions_df(*, predictions_csv_path: str, obs_names: list[str]) -> pd.DataFrame:
    """Load predictions, validate barcode coverage, and order rows by observations."""
    predictions_df = pd.read_csv(predictions_csv_path, index_col=0)
    predictions_df.index = predictions_df.index.astype(str)
    missing_barcodes = set(obs_names) - set(predictions_df.index)
    if missing_barcodes:
        sample = sorted(missing_barcodes)[:5]
        raise ValueError(
            f"{len(missing_barcodes)} barcodes from h5ad not found in predictions CSV "
            f"(first {len(sample)}): {sample}"
        )
    return predictions_df.reindex(obs_names)


def _validate_predictions(
    *,
    predictions_df: pd.DataFrame,
    top_k: int,
    crosswalk_map: dict[str, str],
    known_cl_ids: set[str],
    comparator_name: str,
) -> None:
    """Reject labels missing from crosswalk or ontology."""
    unmapped: set[str] = set()
    non_ontology: set[str] = set()
    for rank in range(1, top_k + 1):
        label_col = f"predicted_label_{rank}"
        for label in predictions_df[label_col].dropna().astype(str).unique():
            cl_id = crosswalk_map.get(label)
            if cl_id is None:
                unmapped.add(label)
            elif cl_id not in known_cl_ids:
                non_ontology.add(cl_id)
    if unmapped:
        sample = sorted(unmapped)[:20]
        raise ValueError(
            f"{comparator_name} class labels missing from the crosswalk: {sample} ({len(unmapped)}); "
            f"CL ids absent from the ontology resource: {sorted(non_ontology)[:20]} ({len(non_ontology)})"
        )
    if non_ontology:
        sample = sorted(non_ontology)[:20]
        raise ValueError(
            f"{comparator_name} class labels missing from the crosswalk: [] (0); "
            f"CL ids absent from the ontology resource: {sample} ({len(non_ontology)})"
        )


def _build_cell_row(
    *,
    barcode: str,
    predictions_df: pd.DataFrame,
    top_k: int,
    crosswalk_map: dict[str, str],
    crosswalk_label_map: dict[str, str],
) -> dict[str, Any]:
    """Build one CAS row, merging duplicate Cell Ontology IDs."""
    entries: dict[str, list[Any]] = {}
    for rank in range(1, top_k + 1):
        label = predictions_df.at[barcode, f"predicted_label_{rank}"]
        if label is None or pd.isna(label):
            continue
        label_str = str(label)
        cl_id = crosswalk_map[label_str]
        cl_label = crosswalk_label_map.get(label_str, cl_id)
        score = predictions_df.at[barcode, f"predicted_score_{rank}"]
        score = float(score) if score is not None and not pd.isna(score) else 0.0
        if cl_id in entries:
            entries[cl_id][1] += score
        else:
            entries[cl_id] = [cl_label, score]

    merged = sorted(entries.items(), key=lambda item: item[1][1], reverse=True)
    row: dict[str, Any] = {}
    for rank, (cl_id, (cl_label, score)) in enumerate(merged[:top_k], start=1):
        row[f"cas_cell_type_label_{rank}"] = cl_label
        row[f"cas_cell_type_name_{rank}"] = cl_id
        row[f"cas_cell_type_score_{rank}"] = score
    return row


def map_ranked_predictions_to_cas_labels(
    *,
    predictions_csv_path: str,
    h5ad_path: str,
    output_dir: str,
    crosswalk_csv_path: str,
    crosswalk_label_col: str,
    crosswalk_cl_id_col: str,
    ontology_resource_path: str,
    model_name: str,
    comparator_name: str,
    crosswalk_cl_label_col: str | None = None,
    top_k: int = 5,
) -> dict[str, str]:
    """Convert ranked predictions into CAS labels and metadata.

    :param predictions_csv_path: Read ranked prediction CSV from this path.
    :param h5ad_path: Read source observations and their order from this H5AD.
    :param output_dir: Write CAS labels and metadata into this directory.
    :param crosswalk_csv_path: Read comparator-label to Cell Ontology mappings from this CSV.
    :param crosswalk_label_col: Read source labels from this crosswalk column.
    :param crosswalk_cl_id_col: Read Cell Ontology IDs from this crosswalk column.
    :param ontology_resource_path: Validate Cell Ontology IDs against this resource.
    :param model_name: Write this benchmark group name into metadata.
    :param comparator_name: Use this comparator name in validation errors.
    :param crosswalk_cl_label_col: Read display labels from this optional crosswalk column.
    :param top_k: Read and write this many prediction ranks.

    :raises ValueError: If barcodes, labels, or Cell Ontology IDs are invalid.
    :returns: Return paths to inferred labels and metadata.
    """
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    adata = anndata.read_h5ad(h5ad_path)
    obs_names = adata.obs_names.astype(str).tolist()
    predictions_df = _load_and_align_predictions_df(
        predictions_csv_path=predictions_csv_path,
        obs_names=obs_names,
    )
    crosswalk_map, crosswalk_label_map = load_crosswalk_maps(
        crosswalk_csv_path=crosswalk_csv_path,
        crosswalk_label_col=crosswalk_label_col,
        crosswalk_cl_id_col=crosswalk_cl_id_col,
        crosswalk_cl_label_col=crosswalk_cl_label_col,
    )

    with open(ontology_resource_path) as file:
        ontology_resource = json.load(file)
    cell_ontology_cache = CellOntologyCache(ontology_resource)
    _validate_predictions(
        predictions_df=predictions_df,
        top_k=top_k,
        crosswalk_map=crosswalk_map,
        known_cl_ids=set(cell_ontology_cache.cl_names_to_idx_map),
        comparator_name=comparator_name,
    )

    label_cols = [f"cas_cell_type_label_{rank}" for rank in range(1, top_k + 1)]
    name_cols = [f"cas_cell_type_name_{rank}" for rank in range(1, top_k + 1)]
    score_cols = [f"cas_cell_type_score_{rank}" for rank in range(1, top_k + 1)]
    rows = [
        _build_cell_row(
            barcode=barcode,
            predictions_df=predictions_df,
            top_k=top_k,
            crosswalk_map=crosswalk_map,
            crosswalk_label_map=crosswalk_label_map,
        )
        for barcode in obs_names
    ]
    labels_df = pd.DataFrame(rows, index=pd.Index(obs_names, name=None), columns=label_cols + name_cols + score_cols)
    labels_path = output_dir_path / "inferred_labels.csv"
    labels_df.to_csv(labels_path)

    metadata = {
        "input_path": str(Path(h5ad_path).resolve()),
        "model_name": model_name,
        "n_cells": len(obs_names),
        "ontology_resource_name": None,
    }
    metadata_path = output_dir_path / "metadata.json"
    with open(metadata_path, "w") as file:
        json.dump(metadata, file, indent=2)

    return {
        "inferred_labels_path": str(labels_path),
        "metadata_path": str(metadata_path),
    }
