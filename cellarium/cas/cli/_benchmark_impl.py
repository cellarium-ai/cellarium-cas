"""
Implementation helpers for ``cellarium-cas benchmark`` commands.

This module imports optional benchmark dependencies at module import time.  Keep
``cellarium.cas.cli.benchmark`` lightweight so Sphinx can import Click command
objects without requiring ``cellarium-cas[benchmark]``.

Pipeline steps
--------------
1. ``run_confusion_matrix_step`` -- build per-sample confusion matrices -> ``cm_raw_k{n}/``
2. ``run_aggregate_step``         -- aggregate per-model               -> ``cm_aggregate_k{n}/``
3. ``run_f_measure_step``         -- compute flat F-measure CSVs (wide format, one column set per k)
4. ``run_hierarchical_f_measure_step`` -- compute hierarchical F-measure CSVs (always uses k=1)
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
    # ONTOLOGY_RESPONSE_FILENAME is checked for existence only (as a directory-validity sentinel);
    # its content is not read by the benchmark pipeline.
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


def _load_sample_data(
    annotate_dir: Path,
    gt_label: str,
    inferred_label: str,
    metadata: t.Dict[str, t.Any],
) -> t.Tuple[t.List[str], t.List[str], pd.DataFrame, t.List[str]]:
    """Load ground-truth labels, top-1 predictions, the full inferred-labels DataFrame, and the label order."""
    input_path = metadata["input_path"]

    adata = anndata.read_h5ad(input_path)
    if gt_label not in adata.obs.columns:
        raise ValueError(
            f"Ground-truth column '{gt_label}' not found in {input_path}.\n"
            f"Available obs columns: {sorted(adata.obs.columns.tolist())}"
        )
    y_true = adata.obs[gt_label].astype(str).tolist()

    labels_csv = annotate_dir / INFERRED_LABELS_FILENAME
    labels_df = pd.read_csv(labels_csv, index_col=0)
    if inferred_label not in labels_df.columns:
        raise ValueError(
            f"Inferred-label column '{inferred_label}' not found in {labels_csv}.\n"
            f"Available columns: {sorted(labels_df.columns.tolist())}"
        )
    y_pred = labels_df[inferred_label].astype(str).tolist()

    ontology_resource = load_ontology_resource(annotate_dir)
    label_order: t.List[str] = ontology_resource["cl_names"]

    return y_true, y_pred, labels_df, label_order


def _build_f_measure_predictions(
    y_true: t.List[str],
    y_pred: t.List[str],
    labels_df: pd.DataFrame,
    inferred_labels_for_k: t.List[str],
) -> t.List[str]:
    """Return per-cell predictions for flat F-measure given a list of ranked label columns.

    When only one column is provided (k=1) the top-1 predictions are returned unchanged.
    Otherwise, a cell's true label replaces its prediction if it appears anywhere in the
    top-k columns; the top-1 prediction is used as fallback.
    """
    if len(inferred_labels_for_k) == 1:
        return y_pred
    ranked_predictions = labels_df[inferred_labels_for_k].astype(str)
    return [
        true_label if true_label in top_k_predictions else fallback_prediction
        for true_label, top_k_predictions, fallback_prediction in zip(
            y_true,
            ranked_predictions.itertuples(index=False, name=None),
            y_pred,
        )
    ]


def _effective_inferred_labels(
    inferred_label: str,
    k: int,
    labels_df: pd.DataFrame,
) -> t.List[str]:
    """Return ranked label columns for k, clamped to those that actually exist in labels_df.

    Stops at the first missing rank so that models with fewer ranked columns than k are
    handled gracefully: their metrics for k >= max_available_rank equal those at max_available_rank.
    """
    all_cols = _ranked_inferred_label_columns(inferred_label, k)
    effective = []
    for col in all_cols:
        if col in labels_df.columns:
            effective.append(col)
        else:
            break
    if len(effective) < len(all_cols):
        logger.debug(
            "k=%d: clamped to %d ranked column(s) — '%s' not found in inferred_labels.csv",
            k,
            len(effective),
            all_cols[len(effective)],
        )
    return effective


def _process_annotate_dir(
    annotate_dir: Path,
    k_cm_dirs: t.Dict[int, Path],
    y_true: t.List[str],
    y_pred: t.List[str],
    labels_df: pd.DataFrame,
    label_order: t.List[str],
    meta_base: t.Dict[str, t.Any],
    inferred_label: str,
) -> None:
    """Build and save one confusion matrix per k value for a single annotate directory."""
    for k, cm_dir in k_cm_dirs.items():
        effective_cols = _effective_inferred_labels(inferred_label, k, labels_df)
        f_measure_y_pred = _build_f_measure_predictions(y_true, y_pred, labels_df, effective_cols)
        cm = build_confusion_matrix(y_true, f_measure_y_pred, label_order)
        meta: t.Dict[str, t.Any] = {
            **meta_base,
            "matrix_shape": list(cm.shape),
            "f_measure_top_k": k,
            "f_measure_inferred_labels": effective_cols,
        }
        sample_cm_dir = cm_dir / annotate_dir.name
        save_confusion_matrix(cm, meta, sample_cm_dir)
        logger.info("Saved k=%d confusion matrix (%d x %d) -> %s", k, cm.shape[0], cm.shape[1], sample_cm_dir)


def run_confusion_matrix_step(
    annotate_dirs: t.Union[str, Path],
    output_dir: t.Union[str, Path],
    gt_label: str,
    inferred_label: str,
    f_measure_top_k: int = 1,
) -> t.Dict[str, t.Any]:
    """
    Build per-sample confusion matrices and write them to ``<output_dir>/cm_raw_k{n}/``.

    One directory is created for each k in 1..f_measure_top_k.  ``cm_raw_k1/`` contains
    top-1 predictions and is the source of truth for hierarchical F-measure.

    :param annotate_dirs: Path to a parent directory whose immediate subdirectories
        are annotate output dirs, or a ``.txt`` file listing one dir per line.
    :param output_dir: Benchmarking workspace directory (will be created if absent).
    :param gt_label: Column name in the original ``.h5ad`` ``obs`` that contains
        ground-truth cell type labels.
    :param inferred_label: Column name in ``inferred_labels.csv`` that contains the
        top-ranked predicted cell type labels (e.g. ``cas_cell_type_name_1``).
    :param f_measure_top_k: Maximum k to evaluate.  Metrics are computed for every k
        from 1 to this value.  Hierarchical F-measure always uses k=1.
    :returns: Dict with ``n_samples`` (int) and ``cm_raw_k1_dir`` (str).
    """
    output_dir_path = Path(output_dir).resolve()

    k_cm_dirs: t.Dict[int, Path] = {}
    for k in range(1, f_measure_top_k + 1):
        d = output_dir_path / f"cm_raw_k{k}"
        d.mkdir(parents=True, exist_ok=True)
        k_cm_dirs[k] = d

    dirs = collect_annotate_output_dirs(
        annotate_dirs,
        required_files=_BENCHMARK_REQUIRED_FILES,
        missing_hint=("Ensure 'cellarium-cas annotate' was run with --infer-labels and --save-ontology-resource."),
    )
    logger.info("Building confusion matrices for %d annotate output directories.", len(dirs))

    for annotate_dir in dirs:
        metadata = load_metadata(annotate_dir)
        model_name = metadata.get("model_name", "unknown")
        logger.info("Processing: %s  (model=%s)", annotate_dir, model_name)

        y_true, y_pred, labels_df, label_order = _load_sample_data(annotate_dir, gt_label, inferred_label, metadata)

        meta_base: t.Dict[str, t.Any] = {
            "model_name": model_name,
            "test_sample": metadata["input_path"],
            "annotate_dir": str(annotate_dir.resolve()),
            "gt_label": gt_label,
            "inferred_label": inferred_label,
            "label_order": label_order,
        }
        _process_annotate_dir(
            annotate_dir, k_cm_dirs, y_true, y_pred, labels_df, label_order, meta_base, inferred_label
        )

    return {"n_samples": len(dirs), "cm_raw_k1_dir": str(k_cm_dirs[1])}


def run_aggregate_step(output_dir: t.Union[str, Path]) -> t.Dict[str, t.Any]:
    """
    Aggregate per-sample confusion matrices by model name.

    Discovers all ``cm_raw_k{n}/`` directories in *output_dir*, groups their sub-directories
    by ``model_name``, sums matrices within each group, and writes aggregated matrices to
    ``cm_aggregate_k{n}/``.

    :param output_dir: Benchmarking workspace directory (must already contain ``cm_raw_k1/``).
    :returns: Dict with ``n_groups`` (int) and ``cm_aggregate_k1_dir`` (str).
    """
    output_dir_path = Path(output_dir).resolve()
    k_raw_dirs = _discover_k_dirs(output_dir_path, "cm_raw_k")

    if not k_raw_dirs:
        raise FileNotFoundError(
            f"No cm_raw_k{{n}} directories found in {output_dir_path}. Run the confusion-matrix step first."
        )

    n_groups = 0
    for k, cm_raw_k_dir in sorted(k_raw_dirs.items()):
        cm_aggregate_k_dir = output_dir_path / f"cm_aggregate_k{k}"
        n_groups = _aggregate_cm_tree(cm_raw_k_dir, cm_aggregate_k_dir)
        logger.info("Saved aggregate CMs for k=%d -> %s", k, cm_aggregate_k_dir)

    return {"n_groups": n_groups, "cm_aggregate_k1_dir": str(output_dir_path / "cm_aggregate_k1")}


def run_f_measure_step(output_dir: t.Union[str, Path]) -> t.Dict[str, t.Any]:
    """
    Compute standard F-measure metrics from per-sample and aggregated confusion matrices.

    Reads all ``cm_raw_k{n}/`` and ``cm_aggregate_k{n}/`` directories and writes a single
    wide-format CSV per granularity level.  Each metric column is suffixed with ``_k{n}``.

    - ``f_measure_per_sample.csv``
    - ``f_measure_per_group.csv``

    :param output_dir: Benchmarking workspace directory.
    :returns: Dict with ``per_sample_path`` and ``per_group_path``.
    """
    output_dir_path = Path(output_dir).resolve()
    k_raw_dirs = _discover_k_dirs(output_dir_path, "cm_raw_k")
    k_aggregate_dirs = _discover_k_dirs(output_dir_path, "cm_aggregate_k")

    if not k_raw_dirs or not k_aggregate_dirs:
        raise FileNotFoundError(
            f"cm_raw_k{{n}} or cm_aggregate_k{{n}} directories not found in {output_dir_path}. "
            "Run 'benchmark confusion-matrix' and 'benchmark aggregate' first."
        )

    # --- per sample ---
    per_sample_rows: t.List[t.Dict[str, t.Any]] = []
    for sample_dir in _iter_cm_dirs(k_raw_dirs[1]):
        _, meta_k1 = load_confusion_matrix(sample_dir)
        row: t.Dict[str, t.Any] = {"model_name": meta_k1["model_name"], "test_sample": meta_k1["test_sample"]}
        for k, raw_k_dir in sorted(k_raw_dirs.items()):
            cm, _ = load_confusion_matrix(raw_k_dir / sample_dir.name)
            metrics = compute_f_measure_from_cm(cm)
            row.update({f"{key}_k{k}": val for key, val in metrics.items()})
        per_sample_rows.append(row)

    per_sample_path = output_dir_path / "f_measure_per_sample.csv"
    pd.DataFrame(per_sample_rows).to_csv(per_sample_path, index=False)
    logger.info("Wrote %s (%d rows)", per_sample_path, len(per_sample_rows))

    # --- per group ---
    per_group_rows: t.List[t.Dict[str, t.Any]] = []
    for group_dir in _iter_cm_dirs(k_aggregate_dirs[1]):
        _, meta_k1 = load_confusion_matrix(group_dir)
        row = {"group_name": meta_k1["model_name"]}
        for k, agg_k_dir in sorted(k_aggregate_dirs.items()):
            cm, _ = load_confusion_matrix(agg_k_dir / group_dir.name)
            metrics = compute_f_measure_from_cm(cm)
            row.update({f"{key}_k{k}": val for key, val in metrics.items()})
        per_group_rows.append(row)

    per_group_path = output_dir_path / "f_measure_per_group.csv"
    pd.DataFrame(per_group_rows).to_csv(per_group_path, index=False)
    logger.info("Wrote %s (%d rows)", per_group_path, len(per_group_rows))

    return {"per_sample_path": str(per_sample_path), "per_group_path": str(per_group_path)}


def run_hierarchical_f_measure_step(output_dir: t.Union[str, Path]) -> t.Dict[str, t.Any]:
    """
    Compute hierarchical F-measure metrics from per-sample and aggregated confusion matrices.

    Always uses ``cm_raw_k1/`` and ``cm_aggregate_k1/`` (top-1 predictions) regardless of
    the ``--f-measure-top-k`` value used during the confusion-matrix step.  Writes:

    - ``hierarchical_f_measure_per_sample.csv``
    - ``hierarchical_f_measure_per_group.csv``

    :param output_dir: Benchmarking workspace directory.
    :returns: Dict with ``per_sample_path`` and ``per_group_path``.
    """
    output_dir_path = Path(output_dir).resolve()
    cm_raw_k1 = output_dir_path / "cm_raw_k1"
    cm_aggregate_k1 = output_dir_path / "cm_aggregate_k1"
    for d in (cm_raw_k1, cm_aggregate_k1):
        if not d.exists():
            raise FileNotFoundError(
                f"Required benchmark directory not found: {d}.\n"
                "Run 'benchmark confusion-matrix' and 'benchmark aggregate' first."
            )

    ontology_cache = _load_ontology_cache(output_dir_path)

    # --- per sample ---
    per_sample_rows: t.List[t.Dict[str, t.Any]] = []
    for sample_dir in _iter_cm_dirs(cm_raw_k1):
        cm, meta = load_confusion_matrix(sample_dir)
        label_order: t.List[str] = meta["label_order"]
        metrics = compute_hierarchical_f_measure_from_cm(cm, label_order, ontology_cache)
        per_sample_rows.append({"model_name": meta["model_name"], "test_sample": meta["test_sample"], **metrics})

    per_sample_path = output_dir_path / "hierarchical_f_measure_per_sample.csv"
    pd.DataFrame(per_sample_rows).to_csv(per_sample_path, index=False)
    logger.info("Wrote %s (%d rows)", per_sample_path, len(per_sample_rows))

    # --- per group ---
    per_group_rows: t.List[t.Dict[str, t.Any]] = []
    for group_dir in _iter_cm_dirs(cm_aggregate_k1):
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
    :param f_measure_top_k: Maximum k for flat F-measure.  Hierarchical F-measure always uses k=1.
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


def _discover_k_dirs(output_dir: Path, prefix: str) -> t.Dict[int, Path]:
    """Return a sorted dict mapping k -> Path for directories named {prefix}{k}/ in output_dir."""
    result: t.Dict[int, Path] = {}
    if not output_dir.exists():
        return result
    for p in output_dir.iterdir():
        if p.is_dir() and p.name.startswith(prefix):
            try:
                k = int(p.name[len(prefix) :])
                result[k] = p
            except ValueError:
                continue
    return dict(sorted(result.items()))


def _aggregate_cm_tree(cm_raw_dir: Path, cm_aggregate_dir: Path) -> int:
    cm_aggregate_dir.mkdir(parents=True, exist_ok=True)

    if not cm_raw_dir.exists():
        raise FileNotFoundError(f"cm_raw directory not found: {cm_raw_dir}. Run the confusion-matrix step first.")

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


def _load_ontology_cache(output_dir: Path) -> CellOntologyCache:
    """Load CellOntologyCache from the first source annotate dir referenced in cm_raw_k1 metadata."""
    cm_raw_k1 = output_dir / "cm_raw_k1"
    for sample_dir in _iter_cm_dirs(cm_raw_k1):
        _, meta = load_confusion_matrix(sample_dir)
        annotate_dir = Path(meta["annotate_dir"])
        resource = load_ontology_resource(annotate_dir)
        return CellOntologyCache(resource)
    raise ValueError(f"No confusion matrix artifacts found in {cm_raw_k1}. Run the confusion-matrix step first.")
