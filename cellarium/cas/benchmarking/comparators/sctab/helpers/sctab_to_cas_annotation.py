"""
Produce a complete CAS-compatible annotate output directory from scTab top-k predictions
in a single call.

Combines :func:`map_sctab_to_cas_labels` and :func:`build_ontology_response` into one
pipeline that writes the same four files produced by ``cellarium-cas annotate``:

- ``inferred_labels.csv``
- ``metadata.json``
- ``ontology_response.json``
- ``ontology_resource.json``
"""

import typing as t

from cellarium.cas.benchmarking.comparators.shared.cas_labels import build_ontology_response

from .map_sctab_to_cas_labels import map_sctab_to_cas_labels


def sctab_to_cas_annotation(
    predictions_csv_path: str,
    h5ad_path: str,
    crosswalk_csv_path: str,
    crosswalk_label_col: str,
    crosswalk_cl_id_col: str,
    ontology_resource_path: str,
    model_name: str,
    output_dir: str,
    crosswalk_cl_label_col: t.Optional[str] = None,
    top_k: int = 5,
) -> t.Dict[str, str]:
    """
    Produce a complete CAS-compatible annotate output directory from scTab predictions.

    Runs :func:`map_sctab_to_cas_labels` to produce ``inferred_labels.csv`` and
    ``metadata.json``, then :func:`build_ontology_response` to produce
    ``ontology_response.json`` and ``ontology_resource.json``.  The resulting *output_dir*
    can be passed directly to ``cellarium-cas benchmark flat`` or
    ``cellarium-cas benchmark ontology-aware``.

    :param predictions_csv_path: Path to the scTab predictions CSV (row index = barcodes,
        columns ``predicted_label_k`` / ``predicted_score_k``).
    :param h5ad_path: Path to the source ``.h5ad`` file; used for authoritative cell ordering.
    :param crosswalk_csv_path: Path to the HRA crosswalk CSV mapping scTab class labels to CL IDs.
    :param crosswalk_label_col: Column in the crosswalk containing scTab class labels.
    :param crosswalk_cl_id_col: Column in the crosswalk containing CL ontology term IDs.
    :param ontology_resource_path: Path to ``ontology_resource.json`` (saved by
        ``cellarium-cas annotate --save-ontology-resource``).  Passed to
        :func:`build_ontology_response`, which loads it and copies it to *output_dir*.
    :param model_name: Benchmark group label (e.g. ``"sctab_tabnet"``).  Written into
        ``metadata.json`` and the ``ontology_response.json`` ``model_name`` field.
    :param output_dir: Directory to write all four output files into (created if absent).
    :param crosswalk_cl_label_col: Optional column in the crosswalk containing human-readable
        CL term labels.  Written to ``cas_cell_type_label_k`` columns.  If ``None``, the CL id
        is used.
    :param top_k: Number of prediction ranks to read from *predictions_csv_path*.

    :returns: Dict with keys ``output_dir``, ``inferred_labels_path``, ``metadata_path``,
        ``ontology_response_path``, and ``ontology_resource_path``.
    """
    from pathlib import Path

    labels_result = map_sctab_to_cas_labels(
        predictions_csv_path=predictions_csv_path,
        h5ad_path=h5ad_path,
        output_dir=output_dir,
        crosswalk_csv_path=crosswalk_csv_path,
        crosswalk_label_col=crosswalk_label_col,
        crosswalk_cl_id_col=crosswalk_cl_id_col,
        ontology_resource_path=ontology_resource_path,
        model_name=model_name,
        crosswalk_cl_label_col=crosswalk_cl_label_col,
        top_k=top_k,
    )

    build_ontology_response(
        inferred_labels_path=labels_result["inferred_labels_path"],
        ontology_resource_path=ontology_resource_path,
        model_name=model_name,
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
