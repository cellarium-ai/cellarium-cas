"""Build complete CAS-compatible annotations from SCimilarity predictions."""

from __future__ import annotations

from pathlib import Path

from cellarium.cas.benchmarking.comparators.scimilarity.helpers.map_scimilarity_to_cas_labels import (
    map_scimilarity_to_cas_labels,
)
from cellarium.cas.benchmarking.comparators.shared.cas_labels import build_ontology_response


def scimilarity_to_cas_annotation(
    *,
    predictions_csv_path: str,
    h5ad_path: str,
    crosswalk_csv_path: str,
    crosswalk_label_col: str,
    crosswalk_cl_id_col: str,
    ontology_resource_path: str,
    output_dir: str,
    crosswalk_cl_label_col: str | None = None,
    top_k: int = 5,
) -> dict[str, str]:
    """Build four CAS annotation files from SCimilarity predictions.

    :param predictions_csv_path: Read ranked prediction CSV from this path.
    :param h5ad_path: Read source observations and their order from this H5AD.
    :param crosswalk_csv_path: Read SCimilarity-label to Cell Ontology mappings from this CSV.
    :param crosswalk_label_col: Read SCimilarity labels from this crosswalk column.
    :param crosswalk_cl_id_col: Read Cell Ontology IDs from this crosswalk column.
    :param ontology_resource_path: Validate and copy this Cell Ontology resource.
    :param output_dir: Write four CAS annotation files into this directory.
    :param crosswalk_cl_label_col: Read display labels from this optional crosswalk column.
    :param top_k: Read and write this many prediction ranks.

    :raises ValueError: If source predictions cannot map to valid Cell Ontology IDs.
    :returns: Return paths to all four CAS annotation files.
    """
    labels_result = map_scimilarity_to_cas_labels(
        predictions_csv_path=predictions_csv_path,
        h5ad_path=h5ad_path,
        output_dir=output_dir,
        crosswalk_csv_path=crosswalk_csv_path,
        crosswalk_label_col=crosswalk_label_col,
        crosswalk_cl_id_col=crosswalk_cl_id_col,
        ontology_resource_path=ontology_resource_path,
        crosswalk_cl_label_col=crosswalk_cl_label_col,
        top_k=top_k,
    )
    build_ontology_response(
        inferred_labels_path=labels_result["inferred_labels_path"],
        ontology_resource_path=ontology_resource_path,
        model_name="scimilarity_model_v1.1",
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
