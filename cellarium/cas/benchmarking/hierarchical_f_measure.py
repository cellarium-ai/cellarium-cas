import inspect
import typing as t

import numpy as np
import pandas as pd
from hiclass.metrics import f1 as hiclass_f1
from hiclass.metrics import precision as hiclass_precision
from hiclass.metrics import recall as hiclass_recall

from cellarium.cas.benchmarking.ontology_aware import _precompute_graph
from cellarium.cas.models import CellTypeOntologyAwareResults

_HICLASS_F1_SUPPORTS_ZERO_DIVISION = "zero_division" in inspect.signature(hiclass_f1).parameters

_SUMMARY_COLUMNS = [
    "n_cells",
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
    ground_truth_ancestor_sets: t.List[t.Set[str]], predicted_term_sets: t.List[t.Set[str]]
) -> t.Tuple[float, float, float]:
    # HiClass literature calls predicted term sets "alpha" and ground-truth
    # ancestor sets "beta"; use explicit names in code to avoid swapping them.
    predicted_denominator = sum(len(predicted_term_set) for predicted_term_set in predicted_term_sets)
    ground_truth_denominator = sum(len(ground_truth_set) for ground_truth_set in ground_truth_ancestor_sets)

    y_true = _sets_to_padded_array(ground_truth_ancestor_sets)
    y_pred = _sets_to_padded_array(predicted_term_sets)
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


def _extract_predicted_term_set(annotation: CellTypeOntologyAwareResults.OntologyAwareAnnotation) -> t.Set[str]:
    return {match.cell_type_ontology_term_id for match in annotation.matches if match.cell_type_ontology_term_id != ""}


def _build_class_level_rows(
    ground_truths: t.List[str],
    ground_truth_ancestor_sets: t.List[t.Set[str]],
    predicted_term_sets: t.List[t.Set[str]],
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
            predicted_term_set = predicted_term_sets[i]
            ground_truth_ancestor_set = ground_truth_ancestor_sets[i]
            tp += float(len(predicted_term_set.intersection(ground_truth_ancestor_set)))
            fp += float(len(predicted_term_set.difference(ground_truth_ancestor_set)))
            fn += float(len(ground_truth_ancestor_set.difference(predicted_term_set)))

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
    response: CellTypeOntologyAwareResults,
    ground_truths: t.List[str],
    ontology_resource: t.Dict[str, t.Any],
    class_level: bool = False,
) -> pd.DataFrame:
    """
    Compute hierarchical F-measure metrics from a CAS ontology-aware response.

    HiClass literature calls the per-cell predicted ontology term set ``alpha_i`` and
    the per-cell ground-truth ancestor set ``beta_i``. This implementation uses the
    explicit names ``predicted_term_sets`` and ``ground_truth_ancestor_sets``.

    :param response: CAS ontology-aware response object.
    :param ground_truths: Ground truth CL term IDs positionally aligned with ``response.data``.
    :param ontology_resource: Raw cell ontology resource dict with ``cl_names`` and
        ``children_dictionary`` keys.
    :param class_level: If ``False`` (default), return one summary row. If ``True``, return
        one row per ground-truth class with pooled binary node-count metrics.
    :raises ValueError: If lengths mismatch or any ground truth term is absent from the ontology.
    """
    if len(ground_truths) != len(response.data):
        raise ValueError(
            f"Length mismatch: ground_truths has {len(ground_truths)} entries "
            f"but response.data has {len(response.data)} entries."
        )

    all_terms = frozenset(ontology_resource["cl_names"])
    unrecognized = [gt for gt in ground_truths if gt not in all_terms]
    if unrecognized:
        raise ValueError(
            f"The following ground truth terms are not present in the ontology resource: {sorted(set(unrecognized))}. "
            f"Ensure all ground truth labels are valid CL term IDs."
        )

    if not response.data:
        return pd.DataFrame()

    graph = _precompute_graph(ontology_resource)
    ground_truth_ancestor_sets = [set(graph.term_ancestors[gt_term]) for gt_term in ground_truths]
    predicted_term_sets = [_extract_predicted_term_set(annotation) for annotation in response.data]

    class_df = _build_class_level_rows(ground_truths, ground_truth_ancestor_sets, predicted_term_sets)
    if class_level:
        return class_df

    hierarchical_precision, hierarchical_recall, hierarchical_f1 = _compute_micro_hierarchical_metrics(
        ground_truth_ancestor_sets, predicted_term_sets
    )
    summary = {
        "n_cells": len(response.data),
        "micro_hierarchical_precision": hierarchical_precision,
        "micro_hierarchical_recall": hierarchical_recall,
        "micro_hierarchical_f1": hierarchical_f1,
        "macro_hierarchical_precision": float(class_df["hierarchical_precision"].mean()),
        "macro_hierarchical_recall": float(class_df["hierarchical_recall"].mean()),
        "macro_hierarchical_f1": float(class_df["hierarchical_f1"].mean()),
        "macro_weighted_hierarchical_precision": float((class_df["weight"] * class_df["hierarchical_precision"]).sum()),
        "macro_weighted_hierarchical_recall": float((class_df["weight"] * class_df["hierarchical_recall"]).sum()),
        "macro_weighted_hierarchical_f1": float((class_df["weight"] * class_df["hierarchical_f1"]).sum()),
    }
    return pd.DataFrame([summary], columns=_SUMMARY_COLUMNS)
