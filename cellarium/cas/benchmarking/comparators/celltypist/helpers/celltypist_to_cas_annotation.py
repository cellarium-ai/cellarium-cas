"""Build complete CAS-compatible annotations from native CellTypist predictions."""

from __future__ import annotations

import json
from pathlib import Path

import anndata
import numpy as np
import pandas as pd

from cellarium.cas.benchmarking.comparators.celltypist.run import PREDICTIONS_FILENAME, STATS_FILENAME
from cellarium.cas.benchmarking.comparators.shared.cas_labels import build_ontology_response
from cellarium.cas.benchmarking.comparators.shared.ranked_predictions import map_ranked_predictions_to_cas_labels


def _load_stats(*, native_dir: Path, h5ad_path: str, top_k: int) -> dict:
    """Read and validate native stats against the requested source H5AD and rank count."""
    stats_path = native_dir / STATS_FILENAME
    if not stats_path.is_file():
        raise ValueError(f"Missing native artifact: {stats_path}")
    with open(stats_path) as file:
        stats = json.load(file)
    for key in ("input_h5ad", "model_name", "model_sha256", "api_version", "mode", "majority_voting", "top_k"):
        if key not in stats:
            raise ValueError(f"Native stats lack '{key}': {stats_path}")
    if (
        stats["mode"] != "best match"
        or stats["majority_voting"] is not False
        or not stats["model_name"].endswith("_raw")
    ):
        raise ValueError("Only raw per-cell CellTypist output (mode='best match', no majority voting) is supported.")
    if Path(stats["input_h5ad"]).resolve() != Path(h5ad_path).resolve():
        raise ValueError(f"Native output was produced from {stats['input_h5ad']}, not {Path(h5ad_path).resolve()}.")
    if top_k < 1 or top_k > int(stats["top_k"]):
        raise ValueError(f"top_k must be between 1 and the exported {stats['top_k']}; got {top_k}.")
    return stats


def _validate_predictions(*, predictions_path: Path, h5ad_path: str, top_k: int) -> None:
    """Require unique, exactly covering barcodes, all consumed columns, and finite [0, 1] scores."""
    if not predictions_path.is_file():
        raise ValueError(f"Missing native artifact: {predictions_path}")
    predictions = pd.read_csv(predictions_path, index_col=0)
    columns = [f"predicted_{kind}_{rank}" for rank in range(1, top_k + 1) for kind in ("label", "score")]
    missing = [column for column in columns if column not in predictions.columns]
    if missing:
        raise ValueError(f"Native predictions lack columns: {missing}")
    index = pd.Index(predictions.index.astype(str))
    obs_names = anndata.read_h5ad(h5ad_path, backed="r").obs_names.astype(str)
    if not index.is_unique or set(index) != set(obs_names):
        raise ValueError("Native prediction barcodes must be unique and exactly cover the source H5AD.")
    for rank in range(1, top_k + 1):
        scores = predictions.loc[predictions[f"predicted_label_{rank}"].notna(), f"predicted_score_{rank}"]
        values = pd.to_numeric(scores, errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
            raise ValueError(f"predicted_score_{rank} must be finite and within [0, 1] for nonempty labels.")


def _validate_crosswalk(*, path: str, label_col: str, cl_id_col: str, cl_label_col: str | None) -> None:
    """Reject missing columns, blank cells, and conflicting rows for one native label."""
    columns = [label_col, cl_id_col] + ([cl_label_col] if cl_label_col is not None else [])
    crosswalk = pd.read_csv(path, dtype=str, keep_default_na=False)
    for column in columns:
        if column not in crosswalk.columns:
            raise ValueError(
                f"Column '{column}' not found in crosswalk CSV. Available columns: {list(crosswalk.columns)}"
            )
        if (crosswalk[column].str.strip() == "").any():
            raise ValueError(f"Crosswalk column '{column}' contains empty values.")
    conflicts = crosswalk.groupby(label_col)[columns[1:]].nunique().gt(1).any(axis=1)
    if conflicts.any():
        raise ValueError(f"Crosswalk maps one label to multiple targets: {sorted(conflicts[conflicts].index)[:20]}")


def celltypist_to_cas_annotation(
    *,
    native_output_dir: str,
    h5ad_path: str,
    crosswalk_csv_path: str,
    crosswalk_label_col: str,
    crosswalk_cl_id_col: str,
    ontology_resource_path: str,
    output_dir: str,
    crosswalk_cl_label_col: str | None = None,
    top_k: int = 5,
) -> dict[str, str]:
    """Build four CAS annotation files from native CellTypist predictions.

    Labels sharing a Cell Ontology ID among the exported top ``top_k`` merge by maximum score, so
    native rank 1 keeps its score. Scores are sigmoid values, not calibrated probabilities.

    :param native_output_dir: Read ``run_celltypist`` artifacts from this directory.
    :param h5ad_path: Read source observations and their order from this H5AD (must be the one annotated).
    :param crosswalk_csv_path: Read model-specific native-label to Cell Ontology mappings (one ID per label).
    :param crosswalk_label_col: Read native labels from this crosswalk column.
    :param crosswalk_cl_id_col: Read Cell Ontology IDs from this crosswalk column.
    :param ontology_resource_path: Validate and copy this Cell Ontology resource.
    :param output_dir: Write four CAS annotation files into this fresh directory.
    :param crosswalk_cl_label_col: Read display labels from this optional crosswalk column.
    :param top_k: Read and write this many prediction ranks.

    :raises ValueError: If native artifacts, crosswalk, labels, or Cell Ontology IDs are invalid, or
        ``output_dir`` exists and is not empty.
    :returns: Return paths to all four CAS annotation files.
    """
    native_dir = Path(native_output_dir)
    out_dir = Path(output_dir)
    if out_dir.exists() and (not out_dir.is_dir() or any(out_dir.iterdir())):
        raise ValueError(f"Output directory exists and is not empty: {out_dir}")
    stats = _load_stats(native_dir=native_dir, h5ad_path=h5ad_path, top_k=top_k)
    predictions_path = native_dir / PREDICTIONS_FILENAME
    _validate_predictions(predictions_path=predictions_path, h5ad_path=h5ad_path, top_k=top_k)
    _validate_crosswalk(
        path=crosswalk_csv_path,
        label_col=crosswalk_label_col,
        cl_id_col=crosswalk_cl_id_col,
        cl_label_col=crosswalk_cl_label_col,
    )

    labels_result = map_ranked_predictions_to_cas_labels(
        predictions_csv_path=str(predictions_path),
        h5ad_path=h5ad_path,
        output_dir=output_dir,
        crosswalk_csv_path=crosswalk_csv_path,
        crosswalk_label_col=crosswalk_label_col,
        crosswalk_cl_id_col=crosswalk_cl_id_col,
        ontology_resource_path=ontology_resource_path,
        model_name=stats["model_name"],
        comparator_name="CellTypist",
        crosswalk_cl_label_col=crosswalk_cl_label_col,
        top_k=top_k,
        duplicate_score_aggregation="max",
    )
    metadata_path = Path(labels_result["metadata_path"])
    with open(metadata_path) as file:
        metadata = json.load(file)
    metadata.update(
        model_sha256=stats["model_sha256"], api_version=stats["api_version"], duplicate_score_aggregation="max"
    )
    with open(metadata_path, "w") as file:
        json.dump(metadata, file, indent=2)

    build_ontology_response(
        inferred_labels_path=labels_result["inferred_labels_path"],
        ontology_resource_path=ontology_resource_path,
        model_name=stats["model_name"],
        output_dir=output_dir,
    )

    output_dir_path = out_dir.resolve()
    return {
        "output_dir": str(output_dir_path),
        "inferred_labels_path": labels_result["inferred_labels_path"],
        "metadata_path": labels_result["metadata_path"],
        "ontology_response_path": str(output_dir_path / "ontology_response.json"),
        "ontology_resource_path": str(output_dir_path / "ontology_resource.json"),
    }
