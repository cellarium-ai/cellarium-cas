"""Map SCimilarity predictions to CAS-compatible label files."""

from __future__ import annotations

from cellarium.cas.benchmarking.comparators.shared.ranked_predictions import (
    map_ranked_predictions_to_cas_labels,
)


def map_scimilarity_to_cas_labels(
    *,
    predictions_csv_path: str,
    h5ad_path: str,
    output_dir: str,
    crosswalk_csv_path: str,
    crosswalk_label_col: str,
    crosswalk_cl_id_col: str,
    ontology_resource_path: str,
    crosswalk_cl_label_col: str | None = None,
    top_k: int = 5,
) -> dict[str, str]:
    """Convert SCimilarity predictions into CAS labels and metadata.

    :param predictions_csv_path: Read ranked prediction CSV from this path.
    :param h5ad_path: Read source observations and their order from this H5AD.
    :param output_dir: Write CAS labels and metadata into this directory.
    :param crosswalk_csv_path: Read SCimilarity-label to Cell Ontology mappings from this CSV.
    :param crosswalk_label_col: Read SCimilarity labels from this crosswalk column.
    :param crosswalk_cl_id_col: Read Cell Ontology IDs from this crosswalk column.
    :param ontology_resource_path: Validate Cell Ontology IDs against this resource.
    :param crosswalk_cl_label_col: Read display labels from this optional crosswalk column.
    :param top_k: Read and write this many prediction ranks.

    :raises ValueError: If barcodes, labels, or Cell Ontology IDs are invalid.
    :returns: Return paths to inferred labels and metadata.
    """
    return map_ranked_predictions_to_cas_labels(
        predictions_csv_path=predictions_csv_path,
        h5ad_path=h5ad_path,
        output_dir=output_dir,
        crosswalk_csv_path=crosswalk_csv_path,
        crosswalk_label_col=crosswalk_label_col,
        crosswalk_cl_id_col=crosswalk_cl_id_col,
        ontology_resource_path=ontology_resource_path,
        model_name="scimilarity_model_v1.1",
        comparator_name="SCimilarity",
        crosswalk_cl_label_col=crosswalk_cl_label_col,
        top_k=top_k,
    )
