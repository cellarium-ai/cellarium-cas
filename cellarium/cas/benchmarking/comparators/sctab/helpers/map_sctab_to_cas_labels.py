"""Map scTab top-k predictions to CAS-compatible ``inferred_labels.csv`` and ``metadata.json``."""

import json
import typing as t
from pathlib import Path

import anndata
import pandas as pd

from cellarium.cas.benchmarking.comparators.shared.crosswalk import load_crosswalk_maps
from cellarium.cas.postprocessing.cell_ontology.cell_ontology_cache import CellOntologyCache


def _load_and_align_predictions_df(predictions_csv_path: str, obs_names: t.List[str]) -> t.Any:
    """Load the predictions CSV, validate barcode coverage, and reindex to *obs_names* order."""
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
    predictions_df: t.Any,
    top_k: int,
    crosswalk_map: t.Dict[str, str],
    known_cl_ids: t.AbstractSet[str],
) -> None:
    """Raise ``ValueError`` if any predicted label is unmapped or maps outside the ontology."""
    unmapped: t.Set[str] = set()
    non_ontology: t.Set[str] = set()
    for k in range(1, top_k + 1):
        label_col = f"predicted_label_{k}"
        for label in predictions_df[label_col].dropna().astype(str).unique():
            cl_id = crosswalk_map.get(label)
            if cl_id is None:
                unmapped.add(label)
            elif cl_id not in known_cl_ids:
                non_ontology.add(cl_id)
    if unmapped:
        sample = sorted(unmapped)[:20]
        raise ValueError(
            f"scTab class labels missing from the crosswalk: {sample} ({len(unmapped)}); "
            f"CL ids absent from the ontology resource: {sorted(non_ontology)[:20]} ({len(non_ontology)})"
        )
    if non_ontology:
        sample = sorted(non_ontology)[:20]
        raise ValueError(
            f"scTab class labels missing from the crosswalk: [] (0); "
            f"CL ids absent from the ontology resource: {sample} ({len(non_ontology)})"
        )


def _build_cell_row(
    barcode: str,
    predictions_df: t.Any,
    top_k: int,
    crosswalk_map: t.Dict[str, str],
    crosswalk_label_map: t.Dict[str, str],
) -> t.Dict[str, t.Any]:
    """Build one output row dict for *barcode*, merging duplicate CL ids by summing scores."""
    entries: t.Dict[str, t.List[t.Any]] = {}
    for k in range(1, top_k + 1):
        label = predictions_df.at[barcode, f"predicted_label_{k}"]
        if label is None or pd.isna(label):
            continue
        label_str = str(label)
        cl_id = crosswalk_map[label_str]
        cl_label = crosswalk_label_map.get(label_str, cl_id)
        score = predictions_df.at[barcode, f"predicted_score_{k}"]
        score = float(score) if score is not None and not pd.isna(score) else 0.0
        if cl_id in entries:
            entries[cl_id][1] += score
        else:
            entries[cl_id] = [cl_label, score]
    # Rank by merged score descending; stable sort keeps first-seen order on ties.
    merged = sorted(entries.items(), key=lambda item: item[1][1], reverse=True)
    row: t.Dict[str, t.Any] = {}
    for rank, (cl_id, (cl_label, score)) in enumerate(merged[:top_k], start=1):
        row[f"cas_cell_type_label_{rank}"] = cl_label
        row[f"cas_cell_type_name_{rank}"] = cl_id
        row[f"cas_cell_type_score_{rank}"] = score
    return row


def map_sctab_to_cas_labels(
    predictions_csv_path: str,
    h5ad_path: str,
    output_dir: str,
    crosswalk_csv_path: str,
    crosswalk_label_col: str,
    crosswalk_cl_id_col: str,
    ontology_resource_path: str,
    model_name: str,
    crosswalk_cl_label_col: t.Optional[str] = None,
    top_k: int = 5,
) -> t.Dict[str, str]:
    """
    Convert scTab top-k predictions to CAS-compatible ``inferred_labels.csv`` and ``metadata.json``.

    Each cell's ``predicted_label_k`` / ``predicted_score_k`` ranks are mapped to CL
    ontology terms via the crosswalk, and ranks mapping to the same CL id are merged by
    summing their scores.  Output columns follow the CAS annotate convention:
    ``cas_cell_type_label_k`` (human-readable CL label), ``cas_cell_type_name_k`` (CL id),
    ``cas_cell_type_score_k`` (score), one rank per column.

    :param predictions_csv_path: Path to the scTab predictions CSV (row index = barcodes,
        columns ``predicted_label_k`` / ``predicted_score_k``).
    :param h5ad_path: Path to the source ``.h5ad`` file; used for authoritative cell ordering.
    :param output_dir: Directory to write output files into (created if absent).
    :param crosswalk_csv_path: Path to the HRA crosswalk CSV mapping scTab class labels to CL IDs.
    :param crosswalk_label_col: Column in the crosswalk containing scTab class labels.
    :param crosswalk_cl_id_col: Column in the crosswalk containing CL ontology term IDs
        (e.g. ``"CL:0000084"``).  Written to ``cas_cell_type_name_k`` columns.
    :param ontology_resource_path: Path to ``ontology_resource.json`` (saved by
        ``cellarium-cas annotate --save-ontology-resource``); its CL id set is used to
        validate mapped terms.
    :param model_name: Benchmark group label (e.g. ``"sctab_tabnet"``), written to
        ``metadata.json`` and used downstream as the benchmark group key.
    :param crosswalk_cl_label_col: Optional column in the crosswalk containing the
        human-readable CL term label (e.g. ``"cl_label"``).  Written to
        ``cas_cell_type_label_k`` columns.  If ``None``, the CL id is used.
    :param top_k: Number of prediction ranks to read from *predictions_csv_path*.

    :returns: Dict with ``inferred_labels_path`` and ``metadata_path``.
    :raises ValueError: If barcodes in the h5ad are missing from the predictions CSV, if a
        predicted label is absent from the crosswalk, or if a mapped CL id is absent from
        the ontology resource.
    """
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    adata = anndata.read_h5ad(h5ad_path)
    obs_names = adata.obs_names.astype(str).tolist()

    predictions_df = _load_and_align_predictions_df(predictions_csv_path, obs_names)

    crosswalk_map, crosswalk_label_map = load_crosswalk_maps(
        crosswalk_csv_path,
        crosswalk_label_col=crosswalk_label_col,
        crosswalk_cl_id_col=crosswalk_cl_id_col,
        crosswalk_cl_label_col=crosswalk_cl_label_col,
    )

    with open(ontology_resource_path) as f:
        ontology_resource = json.load(f)
    cell_ontology_cache = CellOntologyCache(ontology_resource)
    known_cl_ids = set(cell_ontology_cache.cl_names_to_idx_map.keys())

    _validate_predictions(predictions_df, top_k, crosswalk_map, known_cl_ids)

    label_cols = [f"cas_cell_type_label_{k}" for k in range(1, top_k + 1)]
    name_cols = [f"cas_cell_type_name_{k}" for k in range(1, top_k + 1)]
    score_cols = [f"cas_cell_type_score_{k}" for k in range(1, top_k + 1)]

    rows = [
        _build_cell_row(barcode, predictions_df, top_k, crosswalk_map, crosswalk_label_map) for barcode in obs_names
    ]
    ordered_cols = label_cols + name_cols + score_cols
    labels_df = pd.DataFrame(rows, index=pd.Index(obs_names, name=None), columns=ordered_cols)
    labels_path = output_dir_path / "inferred_labels.csv"
    labels_df.to_csv(labels_path)

    # --- metadata ---
    metadata = {
        "input_path": str(Path(h5ad_path).resolve()),
        "model_name": model_name,
        "n_cells": len(obs_names),
        "ontology_resource_name": None,
    }
    metadata_path = output_dir_path / "metadata.json"
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2)

    return {
        "inferred_labels_path": str(labels_path),
        "metadata_path": str(metadata_path),
    }
