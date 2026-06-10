import typing as t

import pandas as pd
from sklearn.metrics import confusion_matrix, f1_score, precision_score, recall_score

_SUMMARY_COLUMNS = [
    "n_cells",
    "micro_flat_precision",
    "micro_flat_recall",
    "micro_flat_f1",
    "macro_flat_precision",
    "macro_flat_recall",
    "macro_flat_f1",
    "macro_weighted_flat_precision",
    "macro_weighted_flat_recall",
    "macro_weighted_flat_f1",
]

_CLASS_LEVEL_COLUMNS = [
    "ground_truth_class",
    "support",
    "weight",
    "tp",
    "fp",
    "fn",
    "flat_precision",
    "flat_recall",
    "flat_f1",
]


def _safe_divide(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator > 0.0 else 0.0


def extract_predictions_from_adata(
    adata: t.Any,
    column_prefix: str = "cas_cell_type_label",
    top_k: int = 5,
) -> t.List[t.List[str]]:
    """
    Extract ranked cell type predictions from an AnnData object's obs columns.

    Reads columns ``{column_prefix}_1`` through ``{column_prefix}_{top_k}`` from
    ``adata.obs`` and returns them as a list of ranked prediction lists, one per cell.

    :param adata: AnnData object with CAS predictions inserted via postprocessing.
    :param column_prefix: Column name prefix to read from ``adata.obs``
        (e.g. ``"cas_cell_type_label"`` or ``"cas_cell_type_name"``).
    :param top_k: Number of ranked predictions to extract per cell.

    :return: List of length ``n_obs``, each element a list of up to ``top_k`` predicted labels.
    """
    cols = [f"{column_prefix}_{rank}" for rank in range(1, top_k + 1) if f"{column_prefix}_{rank}" in adata.obs.columns]
    return [
        [str(value) for value in row if pd.notna(value)] for row in adata.obs[cols].itertuples(index=False, name=None)
    ]


def _effective_predictions_at_k(
    ground_truths: t.List[str],
    predictions: t.List[t.List[str]],
    k: int,
) -> t.Tuple[t.List[str], t.List[str]]:
    """
    Return (ground_truths, effective_predictions) for top-k evaluation.

    The effective prediction for a cell is the ground truth label if it appears
    anywhere in the top-k predictions (a hit), otherwise the top-1 prediction.
    This is the standard multi-class top-k evaluation convention.

    :param ground_truths: Ground truth label per cell.
    :param predictions: Ranked predictions per cell (outer list = cells).
    :param k: Number of top predictions to consider.

    :return: Tuple of (ground_truth_list, effective_prediction_list).
    """
    effective_preds = []
    for gt, cell_preds in zip(ground_truths, predictions):
        top_k_vals = cell_preds[:k]
        top1 = top_k_vals[0] if top_k_vals else ""
        effective_pred = gt if gt in top_k_vals else top1
        effective_preds.append(effective_pred)
    return ground_truths, effective_preds


def _build_class_level_rows(ground_truths: t.List[str], effective_preds: t.List[str]) -> pd.DataFrame:
    n_cells = len(ground_truths)
    ground_truth_classes = sorted(set(ground_truths))
    labels = sorted(set(ground_truths).union(effective_preds))
    matrix = confusion_matrix(ground_truths, effective_preds, labels=labels)
    label_to_index = {label: i for i, label in enumerate(labels)}

    rows: t.List[t.Dict[str, t.Any]] = []
    for ground_truth_class in ground_truth_classes:
        class_index = label_to_index[ground_truth_class]
        tp = float(matrix[class_index, class_index])
        fp = float(matrix[:, class_index].sum() - tp)
        fn = float(matrix[class_index, :].sum() - tp)
        support = int(matrix[class_index, :].sum())

        precision = _safe_divide(tp, tp + fp)
        recall = _safe_divide(tp, tp + fn)
        f1 = _safe_divide(2.0 * precision * recall, precision + recall)

        rows.append(
            {
                "ground_truth_class": ground_truth_class,
                "support": support,
                "weight": support / n_cells,
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "flat_precision": precision,
                "flat_recall": recall,
                "flat_f1": f1,
            }
        )

    return pd.DataFrame(rows, columns=_CLASS_LEVEL_COLUMNS)


def compute_flat_metrics(
    ground_truths: t.List[str],
    predictions: t.List[t.List[str]],
    top_k: t.Optional[int] = None,
    class_level: bool = False,
) -> pd.DataFrame:
    """
    Compute flat (taxonomy-agnostic) classification metrics at ``top_k``.

    Uses standard multi-class top-k evaluation: a prediction at rank k is considered
    correct if the ground truth label appears anywhere in the top-k predictions.
    The output mirrors ``compute_hierarchical_f_measure_metrics``: summary mode returns
    one row with micro, macro, and support-weighted macro metrics; class-level mode
    returns one row per ground-truth class.

    :param ground_truths: Ground truth label per cell. Must have the same length as ``predictions``.
    :param predictions: Ranked predictions per cell. Each inner list contains labels ordered
        by rank (highest confidence first).
    :param top_k: Evaluate metrics at this k. Defaults to the longest prediction list, so all
        supplied ranks are considered.
    :param class_level: If ``False`` (default), return one summary row. If ``True``, return
        one row per ground-truth class with pooled binary counts.

    :return:
        - Summary (``class_level=False``): DataFrame with columns ``n_cells``,
          ``micro_flat_precision``, ``micro_flat_recall``, ``micro_flat_f1``,
          ``macro_flat_precision``, ``macro_flat_recall``, ``macro_flat_f1``,
          ``macro_weighted_flat_precision``, ``macro_weighted_flat_recall``,
          ``macro_weighted_flat_f1``.
        - Class-level (``class_level=True``): DataFrame with columns
          ``ground_truth_class``, ``support``, ``weight``, ``tp``, ``fp``, ``fn``,
          ``flat_precision``, ``flat_recall``, ``flat_f1``.

    :raises ValueError: If ``len(ground_truths) != len(predictions)``, or if ``top_k`` is not positive.
    """
    if len(ground_truths) != len(predictions):
        raise ValueError(
            f"Length mismatch: ground_truths has {len(ground_truths)} entries "
            f"but predictions has {len(predictions)} entries."
        )

    n_cells = len(ground_truths)
    if n_cells == 0:
        return pd.DataFrame()

    if top_k is None:
        top_k = max(len(p) for p in predictions)
    if top_k < 1:
        raise ValueError("top_k must be positive or inferable from at least one prediction.")

    gts, effective_preds = _effective_predictions_at_k(ground_truths, predictions, top_k)
    class_df = _build_class_level_rows(gts, effective_preds)
    if class_level:
        return class_df

    micro_flat_precision = precision_score(gts, effective_preds, average="micro", zero_division=0)
    micro_flat_recall = recall_score(gts, effective_preds, average="micro", zero_division=0)
    micro_flat_f1 = f1_score(gts, effective_preds, average="micro", zero_division=0)
    summary = {
        "n_cells": n_cells,
        "micro_flat_precision": micro_flat_precision,
        "micro_flat_recall": micro_flat_recall,
        "micro_flat_f1": micro_flat_f1,
        "macro_flat_precision": float(class_df["flat_precision"].mean()),
        "macro_flat_recall": float(class_df["flat_recall"].mean()),
        "macro_flat_f1": float(class_df["flat_f1"].mean()),
        "macro_weighted_flat_precision": float((class_df["weight"] * class_df["flat_precision"]).sum()),
        "macro_weighted_flat_recall": float((class_df["weight"] * class_df["flat_recall"]).sum()),
        "macro_weighted_flat_f1": float((class_df["weight"] * class_df["flat_f1"]).sum()),
    }
    return pd.DataFrame([summary], columns=_SUMMARY_COLUMNS)
