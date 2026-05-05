import typing as t

import pandas as pd
from sklearn.metrics import f1_score, precision_score, recall_score


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
    predictions: t.List[t.List[str]] = []
    for _, row in adata.obs.iterrows():
        cell_preds = []
        for rank in range(1, top_k + 1):
            col = f"{column_prefix}_{rank}"
            if col in adata.obs.columns:
                val = row[col]
                if pd.notna(val):
                    cell_preds.append(str(val))
        predictions.append(cell_preds)
    return predictions


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


def compute_flat_metrics(
    ground_truths: t.List[str],
    predictions: t.List[t.List[str]],
    top_k: t.Optional[int] = None,
    cell_level: bool = False,
) -> pd.DataFrame:
    """
    Compute flat (taxonomy-agnostic) classification metrics at each k from 1 to ``top_k``.

    Uses standard multi-class top-k evaluation: a prediction at rank k is considered
    correct if the ground truth label appears anywhere in the top-k predictions.
    Computes accuracy, macro and weighted F1, precision, and recall via scikit-learn.

    :param ground_truths: Ground truth label per cell. Must have the same length as ``predictions``.
    :param predictions: Ranked predictions per cell. Each inner list contains labels ordered
        by rank (highest confidence first). All inner lists should have at least ``top_k`` entries.
    :param top_k: Evaluate metrics at k=1 through k=``top_k``. Defaults to the length of
        the shortest prediction list.
    :param cell_level: If ``False`` (default), return a summary DataFrame with one row per k
        containing aggregated metrics. If ``True``, return a per-cell DataFrame with
        correct/incorrect flags at each k.

    :return:
        - Summary (``cell_level=False``): DataFrame with a single row and columns
          ``n_cells``, ``top_{k}_accuracy``, ``top_{k}_macro_f1``, ``top_{k}_weighted_f1``,
          ``top_{k}_macro_precision``, ``top_{k}_weighted_precision``, ``top_{k}_macro_recall``,
          ``top_{k}_weighted_recall`` for k in 1..``top_k``.
        - Cell-level (``cell_level=True``): DataFrame with columns
          ``cell_index, ground_truth, k, effective_prediction, correct``.

    :raises ValueError: If ``len(ground_truths) != len(predictions)``.
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
        pred_lengths = [len(p) for p in predictions]
        top_k = min(pred_lengths) if pred_lengths else 1

    if cell_level:
        rows = []
        for k in range(1, top_k + 1):
            _, effective_preds = _effective_predictions_at_k(ground_truths, predictions, k)
            for i, (gt, eff_pred) in enumerate(zip(ground_truths, effective_preds)):
                rows.append(
                    {
                        "cell_index": i,
                        "ground_truth": gt,
                        "k": k,
                        "effective_prediction": eff_pred,
                        "correct": gt == eff_pred,
                    }
                )
        return pd.DataFrame(rows)

    summary: t.Dict[str, t.Any] = {"n_cells": n_cells}
    for k in range(1, top_k + 1):
        gts, effective_preds = _effective_predictions_at_k(ground_truths, predictions, k)
        accuracy = sum(gt in preds[:k] for gt, preds in zip(ground_truths, predictions)) / n_cells
        summary[f"top_{k}_accuracy"] = accuracy
        summary[f"top_{k}_macro_f1"] = f1_score(gts, effective_preds, average="macro", zero_division=0)
        summary[f"top_{k}_weighted_f1"] = f1_score(gts, effective_preds, average="weighted", zero_division=0)
        summary[f"top_{k}_macro_precision"] = precision_score(gts, effective_preds, average="macro", zero_division=0)
        summary[f"top_{k}_weighted_precision"] = precision_score(
            gts, effective_preds, average="weighted", zero_division=0
        )
        summary[f"top_{k}_macro_recall"] = recall_score(gts, effective_preds, average="macro", zero_division=0)
        summary[f"top_{k}_weighted_recall"] = recall_score(
            gts, effective_preds, average="weighted", zero_division=0
        )

    return pd.DataFrame([summary])
