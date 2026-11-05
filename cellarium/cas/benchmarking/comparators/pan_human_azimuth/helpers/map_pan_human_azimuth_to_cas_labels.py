"""Map Pan-Human Azimuth native predictions to CAS-compatible label files."""

from __future__ import annotations

import json
from pathlib import Path

import anndata
import pandas as pd

from cellarium.cas.benchmarking.comparators.shared.crosswalk import load_crosswalk_maps
from cellarium.cas.postprocessing.cell_ontology.cell_ontology_cache import CellOntologyCache

MODEL_NAME = "pan_human_azimuth_v1"
# Rank 1 = most granular, matching the CAS convention and the existing Azimuth comparator.
LEVEL_COLUMNS = ("azimuth_fine", "azimuth_medium", "azimuth_broad")


def map_pan_human_azimuth_to_cas_labels(
    *,
    predictions_csv_path: str,
    h5ad_path: str,
    output_dir: str,
    crosswalk_csv_path: str,
    crosswalk_label_col: str,
    crosswalk_cl_id_col: str,
    ontology_resource_path: str,
    crosswalk_cl_label_col: str | None = None,
) -> dict[str, str]:
    """Convert Pan-Human Azimuth predictions into CAS labels and metadata.

    Rank 1, 2, 3 are ``azimuth_fine``, ``azimuth_medium``, ``azimuth_broad``. Refined levels have
    no native confidence, so ``cas_cell_type_score_k`` is left empty rather than invented.

    :param predictions_csv_path: Read ``pan_human_azimuth_predictions.csv`` from this path.
    :param h5ad_path: Read source observations and their order from this H5AD.
    :param output_dir: Write CAS labels and metadata into this directory.
    :param crosswalk_csv_path: Read Pan-Human Azimuth label to Cell Ontology mappings from this CSV.
    :param crosswalk_label_col: Read Pan-Human Azimuth labels from this crosswalk column.
    :param crosswalk_cl_id_col: Read Cell Ontology IDs from this crosswalk column.
    :param ontology_resource_path: Validate Cell Ontology IDs against this resource.
    :param crosswalk_cl_label_col: Read display labels from this optional crosswalk column.

    :raises ValueError: If barcodes, columns, labels, or Cell Ontology IDs are invalid.
    :returns: Return paths to inferred labels and metadata.
    """
    obs_names = anndata.read_h5ad(h5ad_path, backed="r").obs_names.astype(str).tolist()
    predictions = pd.read_csv(predictions_csv_path, index_col=0)
    predictions.index = predictions.index.astype(str)
    absent = [column for column in LEVEL_COLUMNS if column not in predictions.columns]
    if absent:
        raise ValueError(f"Predictions CSV is missing columns: {absent}.")
    missing_barcodes = set(obs_names) - set(predictions.index)
    if missing_barcodes:
        raise ValueError(
            f"{len(missing_barcodes)} barcodes in the H5AD are missing from the predictions CSV: "
            f"{sorted(missing_barcodes)[:5]}"
        )
    levels = predictions.loc[obs_names, list(LEVEL_COLUMNS)]

    crosswalk_map, crosswalk_label_map = load_crosswalk_maps(
        crosswalk_csv_path=crosswalk_csv_path,
        crosswalk_label_col=crosswalk_label_col,
        crosswalk_cl_id_col=crosswalk_cl_id_col,
        crosswalk_cl_label_col=crosswalk_cl_label_col,
    )
    with open(ontology_resource_path) as file:
        known_cl_ids = set(CellOntologyCache(json.load(file)).cl_names_to_idx_map)

    labels = set(levels.stack().astype(str))
    unmapped = sorted(labels - set(crosswalk_map))
    if unmapped:
        raise ValueError(
            f"Pan-Human Azimuth class labels missing from the crosswalk: {unmapped[:20]} ({len(unmapped)})"
        )
    non_ontology = sorted({crosswalk_map[label] for label in labels} - known_cl_ids)
    if non_ontology:
        raise ValueError(f"CL ids absent from the ontology resource: {non_ontology[:20]} ({len(non_ontology)})")

    out = pd.DataFrame(index=pd.Index(obs_names, name=None))
    for rank, column in enumerate(LEVEL_COLUMNS, start=1):
        values = levels[column].astype(str)
        out[f"cas_cell_type_label_{rank}"] = values.map(lambda v: crosswalk_label_map.get(v, crosswalk_map[v]))
        out[f"cas_cell_type_name_{rank}"] = values.map(crosswalk_map)
        out[f"cas_cell_type_score_{rank}"] = float("nan")
    ranks = range(1, len(LEVEL_COLUMNS) + 1)
    out = out[[f"cas_cell_type_{kind}_{rank}" for kind in ("label", "name", "score") for rank in ranks]]

    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    labels_path = output_dir_path / "inferred_labels.csv"
    out.to_csv(labels_path)
    metadata_path = output_dir_path / "metadata.json"
    with open(metadata_path, "w") as file:
        json.dump(
            {
                "input_path": str(Path(h5ad_path).resolve()),
                "model_name": MODEL_NAME,
                "n_cells": len(obs_names),
                "ontology_resource_name": None,
            },
            file,
            indent=2,
        )
    return {"inferred_labels_path": str(labels_path), "metadata_path": str(metadata_path)}
