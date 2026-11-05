"""Map Pan-Human Azimuth native predictions to CAS-compatible label files via the vendor CL map."""

from __future__ import annotations

import json
from importlib.util import find_spec
from pathlib import Path

import anndata
import pandas as pd

from cellarium.cas.benchmarking.comparators.shared.crosswalk import load_crosswalk_maps
from cellarium.cas.postprocessing.cell_ontology.cell_ontology_cache import CellOntologyCache

MODEL_NAME = "pan_human_azimuth_v1"
# Rank 1 = most granular, matching the CAS convention and the existing Azimuth comparator.
LEVEL_COLUMNS = ("azimuth_fine", "azimuth_medium", "azimuth_broad")
VENDOR_LABEL_COL, VENDOR_CL_ID_COL, VENDOR_CL_LABEL_COL = "Annotation_Label", "CL_ID", "CL_Label"


def _vendor_map_path() -> Path:
    """Locate the vendor v1 ``cell_ontology_map.csv`` without importing ``panhumanpy`` (TensorFlow)."""
    spec = find_spec("panhumanpy")
    if spec is None or not spec.submodule_search_locations:
        raise ValueError("panhumanpy is not installed; pass vendor_map_csv_path.")
    return Path(spec.submodule_search_locations[0]) / "_tools" / "v1" / "cell_ontology_map" / "cell_ontology_map.csv"


def map_pan_human_azimuth_to_cas_labels(
    *,
    predictions_csv_path: str,
    h5ad_path: str,
    output_dir: str,
    ontology_resource_path: str,
    vendor_map_csv_path: str | None = None,
    crosswalk_csv_path: str | None = None,
    crosswalk_source_cl_id_col: str | None = None,
    crosswalk_target_cl_id_col: str | None = None,
    crosswalk_cl_label_col: str | None = None,
) -> dict[str, str]:
    """Convert Pan-Human Azimuth predictions into CAS labels and metadata.

    Native labels map to Cell Ontology ids through the vendor's own table. An optional crosswalk then
    translates those CL ids to the benchmark ontology's ids; ids absent from it pass through unchanged.
    Rank 1, 2, 3 are ``azimuth_fine``, ``azimuth_medium``, ``azimuth_broad``. Refined levels have no native
    confidence, so ``cas_cell_type_score_k`` is left empty rather than invented.

    :param predictions_csv_path: Read ``pan_human_azimuth_predictions.csv`` from this path.
    :param h5ad_path: Read source observations and their order from this H5AD.
    :param output_dir: Write CAS labels and metadata into this directory.
    :param ontology_resource_path: Validate final Cell Ontology IDs against this resource.
    :param vendor_map_csv_path: Read the vendor label-to-CL table from this CSV; default is the installed
        ``panhumanpy`` v1 table.
    :param crosswalk_csv_path: Optionally translate vendor CL ids to benchmark CL ids using this CSV.
    :param crosswalk_source_cl_id_col: Read vendor CL ids from this crosswalk column.
    :param crosswalk_target_cl_id_col: Read benchmark CL ids from this crosswalk column.
    :param crosswalk_cl_label_col: Read display labels for the target ids from this optional column.

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
    levels = predictions.loc[obs_names, list(LEVEL_COLUMNS)].astype(str)

    vendor = pd.read_csv(vendor_map_csv_path or _vendor_map_path())
    for col in (VENDOR_LABEL_COL, VENDOR_CL_ID_COL, VENDOR_CL_LABEL_COL):
        if col not in vendor.columns:
            raise ValueError(f"Column '{col}' not found in vendor map. Available columns: {list(vendor.columns)}")
    label_to_cl = dict(zip(vendor[VENDOR_LABEL_COL].astype(str), vendor[VENDOR_CL_ID_COL].astype(str)))
    cl_to_label = dict(zip(vendor[VENDOR_CL_ID_COL].astype(str), vendor[VENDOR_CL_LABEL_COL].astype(str)))

    labels = set(levels.stack())
    unmapped = sorted(labels - set(label_to_cl))
    if unmapped:
        raise ValueError(
            f"Pan-Human Azimuth class labels missing from the vendor CL map: {unmapped[:20]} ({len(unmapped)})"
        )

    if crosswalk_csv_path is not None:
        if crosswalk_source_cl_id_col is None or crosswalk_target_cl_id_col is None:
            raise ValueError("crosswalk_source_cl_id_col and crosswalk_target_cl_id_col are required with a crosswalk.")
        cl_translation, cl_translation_labels = load_crosswalk_maps(
            crosswalk_csv_path=crosswalk_csv_path,
            crosswalk_label_col=crosswalk_source_cl_id_col,
            crosswalk_cl_id_col=crosswalk_target_cl_id_col,
            crosswalk_cl_label_col=crosswalk_cl_label_col,
        )
    else:
        cl_translation, cl_translation_labels = {}, {}

    def final_id(label: str) -> str:
        vendor_id = label_to_cl[label]
        return cl_translation.get(vendor_id, vendor_id)

    def final_label(label: str) -> str:
        vendor_id = label_to_cl[label]
        if vendor_id in cl_translation:
            return cl_translation_labels.get(vendor_id, cl_translation[vendor_id])
        return cl_to_label[vendor_id]

    with open(ontology_resource_path) as file:
        known_cl_ids = set(CellOntologyCache(json.load(file)).cl_names_to_idx_map)
    non_ontology = sorted({final_id(label) for label in labels} - known_cl_ids)
    if non_ontology:
        raise ValueError(
            f"CL ids absent from the ontology resource: {non_ontology[:20]} ({len(non_ontology)}); "
            "add them to the CL crosswalk."
        )

    out = pd.DataFrame(index=pd.Index(obs_names, name=None))
    for rank, column in enumerate(LEVEL_COLUMNS, start=1):
        out[f"cas_cell_type_label_{rank}"] = levels[column].map(final_label)
        out[f"cas_cell_type_name_{rank}"] = levels[column].map(final_id)
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
