"""
Produce a complete CAS-compatible annotate output directory from Azimuth annotations
in a single call.

Combines :func:`map_azimuth_to_cas_labels` and :func:`build_ontology_response` into one
pipeline that writes the same four files produced by ``cellarium-cas annotate``:

- ``inferred_labels.csv``
- ``metadata.json``
- ``ontology_response.json``
- ``ontology_resource.json``
"""

import typing as t

from .build_ontology_response import build_ontology_response
from .map_azimuth_to_cas_labels import map_azimuth_to_cas_labels


def azimuth_to_cas_annotation(
    azimuth_csv_path: str,
    h5ad_path: str,
    crosswalk_csv_path: str,
    crosswalk_azimuth_col: str,
    crosswalk_cl_id_col: str,
    level_specs: t.List[t.Tuple[str, str]],
    azimuth_ref_name: str,
    ontology_resource_path: str,
    output_dir: str,
    crosswalk_cl_name_col: t.Optional[str] = None,
) -> t.Dict[str, str]:
    """
    Produce a complete CAS-compatible annotate output directory from Azimuth annotations.

    Runs :func:`map_azimuth_to_cas_labels` to produce ``inferred_labels.csv`` and
    ``metadata.json``, then :func:`build_ontology_response` to produce
    ``ontology_response.json`` and ``ontology_resource.json``.  The resulting *output_dir*
    can be passed directly to ``cellarium-cas benchmark flat`` or
    ``cellarium-cas benchmark ontology-aware``.

    :param azimuth_csv_path: Path to the Azimuth metadata CSV (row index = barcodes).
    :param h5ad_path: Path to the source ``.h5ad`` file; used for authoritative cell ordering.
    :param crosswalk_csv_path: Path to the HRA crosswalk CSV mapping Azimuth labels to CL IDs.
    :param crosswalk_azimuth_col: Column in the crosswalk containing Azimuth cell type labels.
    :param crosswalk_cl_id_col: Column in the crosswalk containing CL ontology term IDs.
    :param level_specs: List of ``(azimuth_label_col, azimuth_score_col)`` tuples ordered
        **most granular first**.
    :param azimuth_ref_name: Azimuth reference name (e.g. ``"pbmcref"``).  Written into
        ``model_name`` as ``azimuth_<ref_name>``.
    :param ontology_resource_path: Path to ``ontology_resource.json`` (saved by
        ``cellarium-cas annotate --save-ontology-resource``).  Passed to
        :func:`build_ontology_response`, which loads it and copies it to *output_dir*.
    :param output_dir: Directory to write all four output files into (created if absent).
    :param crosswalk_cl_name_col: Optional column in the crosswalk containing human-readable
        CL term names.  Written to ``cas_cell_type_name_k`` columns.  If ``None``, the Azimuth
        label string is used.

    :returns: Dict with keys ``output_dir``, ``inferred_labels_path``, ``metadata_path``,
        ``ontology_response_path``, and ``ontology_resource_path``.
    """
    from pathlib import Path

    labels_result = map_azimuth_to_cas_labels(
        azimuth_csv_path=azimuth_csv_path,
        h5ad_path=h5ad_path,
        output_dir=output_dir,
        crosswalk_csv_path=crosswalk_csv_path,
        crosswalk_azimuth_col=crosswalk_azimuth_col,
        crosswalk_cl_id_col=crosswalk_cl_id_col,
        level_specs=level_specs,
        azimuth_ref_name=azimuth_ref_name,
        crosswalk_cl_name_col=crosswalk_cl_name_col,
    )

    build_ontology_response(
        inferred_labels_path=labels_result["inferred_labels_path"],
        ontology_resource_path=ontology_resource_path,
        azimuth_ref_name=azimuth_ref_name,
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
