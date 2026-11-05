"""Build complete CAS-compatible annotations from Pan-Human Azimuth predictions."""

from __future__ import annotations

from pathlib import Path

from cellarium.cas.benchmarking.comparators.pan_human_azimuth.helpers.map_pan_human_azimuth_to_cas_labels import (
    MODEL_NAME,
    map_pan_human_azimuth_to_cas_labels,
)
from cellarium.cas.benchmarking.comparators.shared.cas_labels import build_ontology_response


def pan_human_azimuth_to_cas_annotation(
    *,
    predictions_csv_path: str,
    h5ad_path: str,
    crosswalk_csv_path: str,
    crosswalk_label_col: str,
    crosswalk_cl_id_col: str,
    ontology_resource_path: str,
    output_dir: str,
    crosswalk_cl_label_col: str | None = None,
) -> dict[str, str]:
    """Build four CAS annotation files from Pan-Human Azimuth predictions.

    :param predictions_csv_path: Read ``pan_human_azimuth_predictions.csv`` from this path.
    :param h5ad_path: Read source observations and their order from this H5AD.
    :param crosswalk_csv_path: Read final-level label to Cell Ontology mappings from this CSV.
    :param crosswalk_label_col: Read Pan-Human Azimuth labels from this crosswalk column.
    :param crosswalk_cl_id_col: Read Cell Ontology IDs from this crosswalk column.
    :param ontology_resource_path: Validate and copy this Cell Ontology resource.
    :param output_dir: Write four CAS annotation files into this directory.
    :param crosswalk_cl_label_col: Read display labels from this optional crosswalk column.

    :raises ValueError: If source predictions cannot map to valid Cell Ontology IDs.
    :returns: Return paths to all four CAS annotation files.
    """
    labels_result = map_pan_human_azimuth_to_cas_labels(
        predictions_csv_path=predictions_csv_path,
        h5ad_path=h5ad_path,
        output_dir=output_dir,
        crosswalk_csv_path=crosswalk_csv_path,
        crosswalk_label_col=crosswalk_label_col,
        crosswalk_cl_id_col=crosswalk_cl_id_col,
        ontology_resource_path=ontology_resource_path,
        crosswalk_cl_label_col=crosswalk_cl_label_col,
    )
    build_ontology_response(
        inferred_labels_path=labels_result["inferred_labels_path"],
        ontology_resource_path=ontology_resource_path,
        model_name=MODEL_NAME,
        output_dir=output_dir,
    )

    output_dir_path = Path(output_dir).resolve()
    return {
        "output_dir": str(output_dir_path),
        "inferred_labels_path": labels_result["inferred_labels_path"],
        "metadata_path": labels_result["metadata_path"],
        "ontology_response_path": str(output_dir_path / "ontology_response.json"),
        "ontology_resource_path": str(output_dir_path / "ontology_resource.json"),
    }
