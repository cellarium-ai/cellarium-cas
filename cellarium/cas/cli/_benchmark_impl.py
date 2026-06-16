"""
Implementation helpers for ``cellarium-cas benchmark`` commands.

This module imports optional benchmark dependencies at module import time.  Keep
``cellarium.cas.cli.benchmark`` lightweight so Sphinx can import Click command
objects without requiring ``cellarium-cas[benchmark]``.

Pipeline steps
--------------
1. ``run_confusion_matrix_step`` -- build per-sample confusion matrices -> ``cm_raw/``
2. ``run_aggregate_step``         -- aggregate per-model               -> ``cm_aggregate/``
3. ``run_f_measure_step``         -- compute F-measure CSVs
4. ``run_hierarchical_f_measure_step`` -- compute hierarchical F-measure CSVs
"""

import typing as t
from pathlib import Path

import anndata
import pandas as pd
import scipy.sparse

from cellarium.cas.benchmarking.confusion_matrix import (
    aggregate_confusion_matrices,
    build_confusion_matrix,
    load_confusion_matrix,
    save_confusion_matrix,
)
from cellarium.cas.benchmarking.f_measure import compute_f_measure_from_cm
from cellarium.cas.benchmarking.hierarchical_f_measure import (
    compute_hierarchical_f_measure_from_cm,
)
from cellarium.cas.logging import logger
from cellarium.cas.postprocessing.cell_ontology.cell_ontology_cache import (
    CellOntologyCache,
)

from ._io import (
    INFERRED_LABELS_FILENAME,
    METADATA_FILENAME,
    ONTOLOGY_RESOURCE_FILENAME,
    ONTOLOGY_RESPONSE_FILENAME,
    collect_annotate_output_dirs,
    load_metadata,
    load_ontology_resource,
)

_BENCHMARK_REQUIRED_FILES = (
    ONTOLOGY_RESPONSE_FILENAME,
    INFERRED_LABELS_FILENAME,
    METADATA_FILENAME,
    ONTOLOGY_RESOURCE_FILENAME,
)


def _ranked_inferred_label_columns(inferred_label: str, top_k: int) -> t.List[str]:
    if top_k == 1:
        return [inferred_label]

    prefix, sep, rank_suffix = inferred_label.rpartition("_")
    if sep == "" or rank_suffix != "1":
        raise ValueError(
            "--f-measure-top-k > 1 requires --inferred-label to end with '_1' so ranked columns can be derived."
        )
    return [f"{prefix}_{rank}" for rank in range(1, top_k + 1)]


def run_confusion_matrix_step(
    annotate_dirs: t.Union[str, Path],
    output_dir: t.Union[str, Path],
    gt_label: str,
    inferred_label: str,
    f_measure_top_k: int = 1,
) -> t.Dict[str, t.Any]:
    """
    Build per-sample confusion matrices and write them to ``<output_dir>/cm_raw/``.

    The top-1 matrices remain the source of truth for hierarchical F-measure.
    Flat F-measure also gets a selected-K matrix tree under
    ``<output_dir>/cm_raw_f_measure/``.

    :param annotate_dirs: Path to a parent directory whose immediate subdirectories
        are annotate output dirs, or a ``.txt`` file listing one dir per line.
    :param output_dir: Benchmarking workspace directory (will be created if absent).
    :param gt_label: Column name in the original ``.h5ad`` ``obs`` that contains
        ground-truth cell type labels.
    :param inferred_label: Column name in ``inferred_labels.csv`` that contains the
        top-ranked predicted cell type labels (e.g. ``cas_cell_type_name_1``).
    :param f_measure_top_k: Number of ranked inferred-label columns to consider
        for flat F-measure.  Hierarchical F-measure remains top-1.
    :returns: Dict with ``n_samples`` (int) and ``cm_raw_dir`` (str).
    """
    output_dir_path = Path(output_dir)
    cm_raw_dir = output_dir_path / "cm_raw"
    cm_raw_f_measure_dir = output_dir_path / "cm_raw_f_measure"
    cm_raw_dir.mkdir(parents=True, exist_ok=True)
    cm_raw_f_measure_dir.mkdir(parents=True, exist_ok=True)

    dirs = collect_annotate_output_dirs(
        annotate_dirs,
        required_files=_BENCHMARK_REQUIRED_FILES,
        missing_hint=("Ensure 'cellarium-cas annotate' was run with --infer-labels and --save-ontology-resource."),
    )
    logger.info("Building confusion matrices for %d annotate output directories.", len(dirs))
    f_measure_inferred_labels = _ranked_inferred_label_columns(inferred_label, f_measure_top_k)

    for annotate_dir in dirs:
        metadata = load_metadata(annotate_dir)
        input_path = metadata["input_path"]
        model_name = metadata.get("model_name", "unknown")
        logger.info("Processing: %s  (model=%s)", annotate_dir, model_name)

        # --- ground truth ---
        adata = anndata.read_h5ad(input_path)
        if gt_label not in adata.obs.columns:
            raise ValueError(
                f"Ground-truth column '{gt_label}' not found in {input_path}.\n"
                f"Available obs columns: {sorted(adata.obs.columns.tolist())}"
            )
        y_true = adata.obs[gt_label].astype(str).tolist()

        # --- predictions ---
        labels_csv = annotate_dir / INFERRED_LABELS_FILENAME
        labels_df = pd.read_csv(labels_csv, index_col=0)
        if inferred_label not in labels_df.columns:
            raise ValueError(
                f"Inferred-label column '{inferred_label}' not found in {labels_csv}.\n"
                f"Available columns: {sorted(labels_df.columns.tolist())}"
            )
        y_pred = labels_df[inferred_label].astype(str).tolist()
        missing_f_measure_columns = [column for column in f_measure_inferred_labels if column not in labels_df.columns]
        if missing_f_measure_columns:
            raise ValueError(
                f"F-measure inferred-label column(s) not found in {labels_csv}: {missing_f_measure_columns}.\n"
                f"Available columns: {labels_df.columns.tolist()}"
            )

        # --- label universe ---
        ontology_resource = load_ontology_resource(annotate_dir)
        label_order: t.List[str] = ontology_resource["cl_names"]

        # --- build + save ---
        cm = build_confusion_matrix(y_true, y_pred, label_order)
        sample_cm_dir = cm_raw_dir / annotate_dir.name
        meta: t.Dict[str, t.Any] = {
            "model_name": model_name,
            "test_sample": input_path,
            "annotate_dir": str(annotate_dir.resolve()),
            "gt_label": gt_label,
            "inferred_label": inferred_label,
            "label_order": label_order,
            "matrix_shape": list(cm.shape),
        }
        save_confusion_matrix(cm, meta, sample_cm_dir)
        logger.info("Saved confusion matrix (%d x %d) -> %s", cm.shape[0], cm.shape[1], sample_cm_dir)

        if f_measure_top_k == 1:
            f_measure_y_pred = y_pred
        else:
            ranked_predictions = labels_df[f_measure_inferred_labels].astype(str)
            f_measure_y_pred = [
                true_label if true_label in top_k_predictions else fallback_prediction
                for true_label, top_k_predictions, fallback_prediction in zip(
                    y_true,
                    ranked_predictions.itertuples(index=False, name=None),
                    y_pred,
                )
            ]

        f_measure_cm = build_confusion_matrix(y_true, f_measure_y_pred, label_order)
        sample_f_measure_cm_dir = cm_raw_f_measure_dir / annotate_dir.name
        f_measure_meta: t.Dict[str, t.Any] = {
            **meta,
            "matrix_shape": list(f_measure_cm.shape),
            "f_measure_top_k": f_measure_top_k,
            "f_measure_inferred_labels": f_measure_inferred_labels,
        }
        save_confusion_matrix(f_measure_cm, f_measure_meta, sample_f_measure_cm_dir)
        logger.info(
            "Saved flat F-measure confusion matrix (%d x %d) -> %s",
            f_measure_cm.shape[0],
            f_measure_cm.shape[1],
            sample_f_measure_cm_dir,
        )

    return {"n_samples": len(dirs), "cm_raw_dir": str(cm_raw_dir)}


def run_aggregate_step(output_dir: t.Union[str, Path]) -> t.Dict[str, t.Any]:
    """
    Aggregate per-sample confusion matrices by model name.

    Reads all sub-directories of ``<output_dir>/cm_raw/``, groups them by the
    ``model_name`` field stored in each matrix's metadata, sums the matrices
    within each group, and writes one aggregated matrix per group to
    ``<output_dir>/cm_aggregate/<model_name>/``.

    :param output_dir: Benchmarking workspace directory (must already contain ``cm_raw/``).
    :returns: Dict with ``n_groups`` (int) and ``cm_aggregate_dir`` (str).
    """
    output_dir_path = Path(output_dir)
    cm_raw_dir = output_dir_path / "cm_raw"
    cm_aggregate_dir = output_dir_path / "cm_aggregate"

    n_groups = _aggregate_cm_tree(cm_raw_dir, cm_aggregate_dir)

    cm_raw_f_measure_dir = output_dir_path / "cm_raw_f_measure"
    if cm_raw_f_measure_dir.exists():
        cm_aggregate_f_measure_dir = output_dir_path / "cm_aggregate_f_measure"
        _aggregate_cm_tree(cm_raw_f_measure_dir, cm_aggregate_f_measure_dir)
        logger.info("Saved flat F-measure aggregate CMs -> %s", cm_aggregate_f_measure_dir)

    return {"n_groups": n_groups, "cm_aggregate_dir": str(cm_aggregate_dir)}


def _aggregate_cm_tree(cm_raw_dir: Path, cm_aggregate_dir: Path) -> int:
    cm_aggregate_dir.mkdir(parents=True, exist_ok=True)

    if not cm_raw_dir.exists():
        raise FileNotFoundError(f"cm_raw directory not found: {cm_raw_dir}. Run the confusion-matrix step first.")

    # Group sample directories by model_name, carrying the matrix shape so we can
    # initialise a zero accumulator without an extra I/O round-trip.
    groups: t.Dict[str, t.List[t.Tuple[Path, t.Tuple[int, int], t.List[str]]]] = {}
    for sample_dir in sorted(p for p in cm_raw_dir.iterdir() if p.is_dir()):
        _, meta = load_confusion_matrix(sample_dir)
        model_name: str = meta["model_name"]
        shape: t.Tuple[int, int] = tuple(meta["matrix_shape"])  # type: ignore[assignment]
        label_order_scan: t.List[str] = meta["label_order"]
        groups.setdefault(model_name, []).append((sample_dir, shape, label_order_scan))

    if not groups:
        raise ValueError(f"No confusion matrix artifacts found in {cm_raw_dir}.")

    logger.info("Aggregating confusion matrices from %s for %d model group(s).", cm_raw_dir, len(groups))

    for model_name, entries in sorted(groups.items()):
        first_shape = entries[0][1]
        label_order: t.List[str] = entries[0][2]

        _validate_group_entries(model_name, entries, first_shape, label_order)

        cms: t.List[scipy.sparse.csr_matrix] = []
        source_annotate_dirs: t.List[str] = []
        n_samples = len(entries)

        for sample_dir, _, _ in entries:
            cm, meta = load_confusion_matrix(sample_dir)
            cms.append(cm)
            source_annotate_dirs.append(meta["annotate_dir"])

        agg_cm = aggregate_confusion_matrices(cms)
        safe_name = _safe_dirname(model_name)
        group_dir = cm_aggregate_dir / safe_name
        agg_meta: t.Dict[str, t.Any] = {
            "model_name": model_name,
            "source_annotate_dirs": source_annotate_dirs,
            "label_order": label_order,
            "matrix_shape": list(agg_cm.shape),
            "n_samples_aggregated": n_samples,
        }
        save_confusion_matrix(agg_cm, agg_meta, group_dir)
        logger.info("Saved aggregated CM for model '%s' (%d samples) -> %s", model_name, n_samples, group_dir)

    return len(groups)


def run_f_measure_step(output_dir: t.Union[str, Path]) -> t.Dict[str, t.Any]:
    """
    Compute standard F-measure metrics from per-sample and aggregated confusion matrices.

    Reads matrices from ``cm_raw/`` and ``cm_aggregate/`` (both must exist) and writes:

    - ``f_measure_per_sample.csv``
    - ``f_measure_per_group.csv``

    :param output_dir: Benchmarking workspace directory.
    :returns: Dict with ``per_sample_path`` and ``per_group_path``.
    """
    output_dir_path = Path(output_dir)
    cm_raw_dir, cm_aggregate_dir = _resolve_f_measure_cm_dirs(output_dir_path)

    # --- per sample ---
    per_sample_rows: t.List[t.Dict[str, t.Any]] = []
    for sample_dir in _iter_cm_dirs(cm_raw_dir):
        cm, meta = load_confusion_matrix(sample_dir)
        metrics = compute_f_measure_from_cm(cm)
        per_sample_rows.append({"model_name": meta["model_name"], "test_sample": meta["test_sample"], **metrics})

    per_sample_path = output_dir_path / "f_measure_per_sample.csv"
    pd.DataFrame(per_sample_rows).to_csv(per_sample_path, index=False)
    logger.info("Wrote %s (%d rows)", per_sample_path, len(per_sample_rows))

    # --- per group ---
    per_group_rows: t.List[t.Dict[str, t.Any]] = []
    for group_dir in _iter_cm_dirs(cm_aggregate_dir):
        cm, meta = load_confusion_matrix(group_dir)
        metrics = compute_f_measure_from_cm(cm)
        per_group_rows.append({"group_name": meta["model_name"], **metrics})

    per_group_path = output_dir_path / "f_measure_per_group.csv"
    pd.DataFrame(per_group_rows).to_csv(per_group_path, index=False)
    logger.info("Wrote %s (%d rows)", per_group_path, len(per_group_rows))

    return {"per_sample_path": str(per_sample_path), "per_group_path": str(per_group_path)}


def run_hierarchical_f_measure_step(output_dir: t.Union[str, Path]) -> t.Dict[str, t.Any]:
    """
    Compute hierarchical F-measure metrics from per-sample and aggregated confusion matrices.

    Uses the ontology ancestor mapping loaded from one of the source annotate directories
    referenced in the per-sample CM metadata.  Writes:

    - ``hierarchical_f_measure_per_sample.csv``
    - ``hierarchical_f_measure_per_group.csv``

    :param output_dir: Benchmarking workspace directory.
    :returns: Dict with ``per_sample_path`` and ``per_group_path``.
    """
    output_dir_path = Path(output_dir)
    _require_cm_dirs(output_dir_path)

    ontology_cache = _load_ontology_cache(output_dir_path)

    # --- per sample ---
    per_sample_rows: t.List[t.Dict[str, t.Any]] = []
    for sample_dir in _iter_cm_dirs(output_dir_path / "cm_raw"):
        cm, meta = load_confusion_matrix(sample_dir)
        label_order: t.List[str] = meta["label_order"]
        metrics = compute_hierarchical_f_measure_from_cm(cm, label_order, ontology_cache)
        per_sample_rows.append({"model_name": meta["model_name"], "test_sample": meta["test_sample"], **metrics})

    per_sample_path = output_dir_path / "hierarchical_f_measure_per_sample.csv"
    pd.DataFrame(per_sample_rows).to_csv(per_sample_path, index=False)
    logger.info("Wrote %s (%d rows)", per_sample_path, len(per_sample_rows))

    # --- per group ---
    per_group_rows: t.List[t.Dict[str, t.Any]] = []
    for group_dir in _iter_cm_dirs(output_dir_path / "cm_aggregate"):
        cm, meta = load_confusion_matrix(group_dir)
        label_order = meta["label_order"]
        metrics = compute_hierarchical_f_measure_from_cm(cm, label_order, ontology_cache)
        per_group_rows.append({"group_name": meta["model_name"], **metrics})

    per_group_path = output_dir_path / "hierarchical_f_measure_per_group.csv"
    pd.DataFrame(per_group_rows).to_csv(per_group_path, index=False)
    logger.info("Wrote %s (%d rows)", per_group_path, len(per_group_rows))

    return {"per_sample_path": str(per_sample_path), "per_group_path": str(per_group_path)}


def run_all_steps(
    annotate_dirs: t.Union[str, Path],
    output_dir: t.Union[str, Path],
    gt_label: str,
    inferred_label: str,
    f_measure_top_k: int = 1,
) -> t.Dict[str, t.Any]:
    """
    Run the full benchmark pipeline: confusion-matrix -> aggregate -> f-measure -> hierarchical.

    :param annotate_dirs: Annotate output directories (parent dir or ``.txt`` list).
    :param output_dir: Benchmarking workspace directory.
    :param gt_label: Ground-truth obs column name in the original ``.h5ad``.
    :param inferred_label: Predicted-label column name in ``inferred_labels.csv``.
    :param f_measure_top_k: Number of ranked inferred-label columns to consider
        for flat F-measure.  Hierarchical F-measure remains top-1.
    :returns: Merged dict of results from all four steps.
    """
    result: t.Dict[str, t.Any] = {}
    result.update(
        run_confusion_matrix_step(
            annotate_dirs,
            output_dir,
            gt_label,
            inferred_label,
            f_measure_top_k=f_measure_top_k,
        )
    )
    result.update(run_aggregate_step(output_dir))
    result.update({"f_measure_" + k: v for k, v in run_f_measure_step(output_dir).items()})
    result.update({"h_f_measure_" + k: v for k, v in run_hierarchical_f_measure_step(output_dir).items()})
    return result


def _validate_group_entries(
    model_name: str,
    entries: t.List[t.Tuple[Path, t.Tuple[int, int], t.List[str]]],
    expected_shape: t.Tuple[int, int],
    expected_label_order: t.List[str],
) -> None:
    """Raise ValueError if any entry's shape or label order differs from the expected values."""
    for _, shape, lo in entries:
        if lo != expected_label_order:
            raise ValueError(
                f"Label order mismatch for model '{model_name}': matrices from different ontology "
                "versions cannot be aggregated. Rerun the confusion-matrix step with a consistent "
                "ontology resource."
            )
        if shape != expected_shape:
            raise ValueError(f"Shape mismatch for model '{model_name}': {shape} vs {expected_shape}.")


def _safe_dirname(name: str) -> str:
    """Replace characters that are problematic on common filesystems."""
    return name.replace("/", "_").replace("\\", "_").replace(":", "_")


def _iter_cm_dirs(parent: Path) -> t.Iterator[Path]:
    """Yield subdirectories of *parent* in sorted order."""
    for p in sorted(parent.iterdir()):
        if p.is_dir():
            yield p


def _resolve_f_measure_cm_dirs(output_dir: Path) -> t.Tuple[Path, Path]:
    cm_raw_f_measure = output_dir / "cm_raw_f_measure"
    cm_aggregate_f_measure = output_dir / "cm_aggregate_f_measure"
    flat_raw_exists = cm_raw_f_measure.exists()
    flat_aggregate_exists = cm_aggregate_f_measure.exists()

    if flat_raw_exists and flat_aggregate_exists:
        return cm_raw_f_measure, cm_aggregate_f_measure

    if flat_raw_exists != flat_aggregate_exists:
        raise FileNotFoundError(
            "Flat F-measure confusion-matrix artifacts are incomplete. Run 'benchmark aggregate' after "
            "'benchmark confusion-matrix'."
        )

    _require_cm_dirs(output_dir)
    return output_dir / "cm_raw", output_dir / "cm_aggregate"


def _require_cm_dirs(output_dir: Path) -> None:
    cm_raw = output_dir / "cm_raw"
    cm_agg = output_dir / "cm_aggregate"
    missing = [str(d) for d in (cm_raw, cm_agg) if not d.exists()]
    if missing:
        raise FileNotFoundError(
            f"Required benchmark subdirectories not found: {missing}.\n"
            "Run 'benchmark confusion-matrix' and 'benchmark aggregate' first."
        )


def _load_ontology_cache(output_dir: Path) -> CellOntologyCache:
    """Load CellOntologyCache from the first source annotate dir referenced in cm_raw metadata."""
    cm_raw = output_dir / "cm_raw"
    for sample_dir in sorted(p for p in cm_raw.iterdir() if p.is_dir()):
        _, meta = load_confusion_matrix(sample_dir)
        annotate_dir = Path(meta["annotate_dir"])
        resource = load_ontology_resource(annotate_dir)
        return CellOntologyCache(resource)
    raise ValueError(f"No confusion matrix artifacts found in {cm_raw}. Run the confusion-matrix step first.")
