import inspect
import typing as t

import numpy as np
import pandas as pd
from hiclass.metrics import f1 as hiclass_f1
from hiclass.metrics import precision as hiclass_precision
from hiclass.metrics import recall as hiclass_recall

from cellarium.cas.benchmarking.ontology_aware import _precompute_graph

_HICLASS_F1_SUPPORTS_ZERO_DIVISION = "zero_division" in inspect.signature(hiclass_f1).parameters

_SUMMARY_METRICS = [
    "micro_hierarchical_precision",
    "micro_hierarchical_recall",
    "micro_hierarchical_f1",
    "macro_hierarchical_precision",
    "macro_hierarchical_recall",
    "macro_hierarchical_f1",
    "macro_weighted_hierarchical_precision",
    "macro_weighted_hierarchical_recall",
    "macro_weighted_hierarchical_f1",
]

_CLASS_LEVEL_COLUMNS = [
    "ground_truth_class",
    "support",
    "weight",
    "tp",
    "fp",
    "fn",
    "hierarchical_precision",
    "hierarchical_recall",
    "hierarchical_f1",
]


def _safe_divide(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator > 0.0 else 0.0


def _sets_to_padded_array(label_sets: t.List[t.Set[str]]) -> np.ndarray:
    max_len = max((len(label_set) for label_set in label_sets), default=0)
    if max_len == 0:
        max_len = 1

    rows = [sorted(label_set) + [""] * (max_len - len(label_set)) for label_set in label_sets]
    return np.asarray(rows, dtype=str)


def _compute_hiclass_f1(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    if _HICLASS_F1_SUPPORTS_ZERO_DIVISION:
        return float(hiclass_f1(y_true, y_pred, average="micro", zero_division=0.0))
    return float(hiclass_f1(y_true, y_pred, average="micro"))


def _compute_micro_hierarchical_metrics(
    ground_truth_ancestor_sets: t.List[t.Set[str]], predicted_ancestor_sets: t.List[t.Set[str]]
) -> t.Tuple[float, float, float]:
    # HiClass literature calls predicted term sets "alpha" and ground-truth
    # ancestor sets "beta"; use explicit names in code to avoid swapping them.
    predicted_denominator = sum(len(predicted_ancestor_set) for predicted_ancestor_set in predicted_ancestor_sets)
    ground_truth_denominator = sum(len(ground_truth_set) for ground_truth_set in ground_truth_ancestor_sets)

    y_true = _sets_to_padded_array(ground_truth_ancestor_sets)
    y_pred = _sets_to_padded_array(predicted_ancestor_sets)
    hierarchical_precision = (
        float(hiclass_precision(y_true, y_pred, average="micro")) if predicted_denominator > 0 else 0.0
    )
    hierarchical_recall = (
        float(hiclass_recall(y_true, y_pred, average="micro")) if ground_truth_denominator > 0 else 0.0
    )
    hierarchical_f1 = (
        _compute_hiclass_f1(y_true, y_pred)
        if predicted_denominator > 0
        and ground_truth_denominator > 0
        and (hierarchical_precision + hierarchical_recall) > 0.0
        else 0.0
    )

    return hierarchical_precision, hierarchical_recall, hierarchical_f1


def _build_predicted_ancestor_set(predicted_terms: t.List[str], graph: t.Any) -> t.Set[str]:
    ancestor_set: t.Set[str] = set()
    for term in predicted_terms:
        if term:
            ancestor_set.update(graph.term_ancestors[term])
    return ancestor_set


def _build_class_level_rows(
    ground_truths: t.List[str],
    ground_truth_ancestor_sets: t.List[t.Set[str]],
    predicted_ancestor_sets: t.List[t.Set[str]],
) -> pd.DataFrame:
    n_cells = len(ground_truths)
    indices_by_class: t.Dict[str, t.List[int]] = {}
    for i, ground_truth_class in enumerate(ground_truths):
        indices_by_class.setdefault(ground_truth_class, []).append(i)

    rows: t.List[t.Dict[str, t.Any]] = []
    for ground_truth_class in sorted(indices_by_class):
        indices = indices_by_class[ground_truth_class]
        tp = 0.0
        fp = 0.0
        fn = 0.0

        for i in indices:
            predicted_ancestor_set = predicted_ancestor_sets[i]
            ground_truth_ancestor_set = ground_truth_ancestor_sets[i]
            tp += float(len(predicted_ancestor_set.intersection(ground_truth_ancestor_set)))
            fp += float(len(predicted_ancestor_set.difference(ground_truth_ancestor_set)))
            fn += float(len(ground_truth_ancestor_set.difference(predicted_ancestor_set)))

        precision = _safe_divide(tp, tp + fp)
        recall = _safe_divide(tp, tp + fn)
        f1 = _safe_divide(2.0 * precision * recall, precision + recall)
        support = len(indices)

        rows.append(
            {
                "ground_truth_class": ground_truth_class,
                "support": support,
                "weight": support / n_cells,
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "hierarchical_precision": precision,
                "hierarchical_recall": recall,
                "hierarchical_f1": f1,
            }
        )

    return pd.DataFrame(rows, columns=_CLASS_LEVEL_COLUMNS)


def compute_hierarchical_f_measure_metrics(
    predictions: t.List[t.List[str]],
    ground_truths: t.List[str],
    ontology_resource: t.Dict[str, t.Any],
    top_k: int = 1,
    class_level: bool = False,
) -> pd.DataFrame:
    """
    Compute hierarchical F-measure metrics from final predicted CL labels.

    Each prediction list contains one or more final predicted CL term IDs for a cell
    (for example, ``cas_cell_type_name_1`` … ``cas_cell_type_name_k`` from
    ``inferred_labels.csv``). Summary mode evaluates every k from 1 to ``top_k`` and
    returns columns prefixed with ``top_{k}_``. Class-level mode returns pooled
    binary node-count metrics for exactly ``top_k`` predictions per cell.

    HiClass literature calls the per-cell predicted ontology term set ``alpha_i`` and
    the per-cell ground-truth ancestor set ``beta_i``. This implementation uses the
    explicit names ``predicted_ancestor_sets`` and ``ground_truth_ancestor_sets``.

    :param predictions: Predicted CL term IDs per cell, positionally aligned with ``ground_truths``.
    :param ground_truths: Ground truth CL term IDs positionally aligned with ``predictions``.
    :param ontology_resource: Raw cell ontology resource dict with ``cl_names`` and
        ``children_dictionary`` keys.
    :param top_k: Highest ranked prediction depth to evaluate. Summary mode returns
        metrics for k=1..``top_k``; class-level mode uses exactly ``top_k``.
    :param class_level: If ``False`` (default), return one summary row. If ``True``, return
        one row per ground-truth class with pooled binary node-count metrics.
    :raises ValueError: If lengths mismatch, ``top_k`` is not positive, or any ground truth/predicted term is absent from the ontology.
    """
    if len(ground_truths) != len(predictions):
        raise ValueError(
            f"Length mismatch: ground_truths has {len(ground_truths)} entries "
            f"but predictions has {len(predictions)} entries."
        )
    if top_k < 1:
        raise ValueError("top_k must be positive.")

    all_terms = frozenset(ontology_resource["cl_names"])
    unrecognized = [gt for gt in ground_truths if gt not in all_terms]
    if unrecognized:
        raise ValueError(
            f"The following ground truth terms are not present in the ontology resource: {sorted(set(unrecognized))}. "
            f"Ensure all ground truth labels are valid CL term IDs."
        )

    unrecognized_predictions = [
        term for cell_predictions in predictions for term in cell_predictions[:top_k] if term and term not in all_terms
    ]
    if unrecognized_predictions:
        raise ValueError(
            "The following predicted terms are not present in the ontology resource: "
            f"{sorted(set(unrecognized_predictions))}. Ensure inferred labels contain valid CL term IDs."
        )

    if not predictions:
        return pd.DataFrame()

    graph = _precompute_graph(ontology_resource)
    ground_truth_ancestor_sets = [set(graph.term_ancestors[gt_term]) for gt_term in ground_truths]

    if class_level:
        predicted_ancestor_sets = [
            _build_predicted_ancestor_set(cell_predictions[:top_k], graph) for cell_predictions in predictions
        ]
        return _build_class_level_rows(ground_truths, ground_truth_ancestor_sets, predicted_ancestor_sets)

    summary: t.Dict[str, t.Any] = {"n_cells": len(predictions)}
    for k in range(1, top_k + 1):
        predicted_ancestor_sets = [
            _build_predicted_ancestor_set(cell_predictions[:k], graph) for cell_predictions in predictions
        ]
        class_df = _build_class_level_rows(ground_truths, ground_truth_ancestor_sets, predicted_ancestor_sets)
        hierarchical_precision, hierarchical_recall, hierarchical_f1 = _compute_micro_hierarchical_metrics(
            ground_truth_ancestor_sets, predicted_ancestor_sets
        )
        summary.update(
            {
                f"top_{k}_micro_hierarchical_precision": hierarchical_precision,
                f"top_{k}_micro_hierarchical_recall": hierarchical_recall,
                f"top_{k}_micro_hierarchical_f1": hierarchical_f1,
                f"top_{k}_macro_hierarchical_precision": float(class_df["hierarchical_precision"].mean()),
                f"top_{k}_macro_hierarchical_recall": float(class_df["hierarchical_recall"].mean()),
                f"top_{k}_macro_hierarchical_f1": float(class_df["hierarchical_f1"].mean()),
                f"top_{k}_macro_weighted_hierarchical_precision": float(
                    (class_df["weight"] * class_df["hierarchical_precision"]).sum()
                ),
                f"top_{k}_macro_weighted_hierarchical_recall": float(
                    (class_df["weight"] * class_df["hierarchical_recall"]).sum()
                ),
                f"top_{k}_macro_weighted_hierarchical_f1": float(
                    (class_df["weight"] * class_df["hierarchical_f1"]).sum()
                ),
            }
        )

    columns = ["n_cells"] + [f"top_{k}_{metric}" for k in range(1, top_k + 1) for metric in _SUMMARY_METRICS]
    return pd.DataFrame([summary], columns=columns)
