"""
Confusion-matrix building, storage, and aggregation for CAS benchmarking.

The confusion matrix is the shared base artifact from which both F-measure and
hierarchical F-measure metrics are derived.  Rows represent true labels; columns
represent predicted labels.  All matrices are expanded to the full Cell Ontology
``cl_names`` label universe so that aggregation across runs is always safe.
"""

import json
import typing as t
from pathlib import Path

import numpy as np
import scipy.sparse
from sklearn.metrics import confusion_matrix as _sklearn_confusion_matrix

CM_MATRIX_FILENAME = "confusion_matrix.npz"
CM_META_FILENAME = "confusion_matrix_meta.json"


def build_confusion_matrix(
    y_true: t.List[str],
    y_pred: t.List[str],
    label_order: t.List[str],
) -> scipy.sparse.csr_matrix:
    """
    Build a sparse confusion matrix aligned to *label_order*.

    The matrix is ``len(label_order) × len(label_order)``.  Row *i* is the true
    class at position *i* in *label_order*; column *j* is the predicted class.

    :param y_true: Ground-truth label per cell.
    :param y_pred: Predicted label per cell (same length as *y_true*).
    :param label_order: Ordered list of all valid labels (the ontology ``cl_names``).
        This determines the row/column ordering of the matrix.
    :raises ValueError: If any label in *y_true* or *y_pred* is not present in
        *label_order*.  Check your crosswalk and ground-truth column.
    :returns: CSR sparse matrix of shape ``(len(label_order), len(label_order))``.
    """
    if len(y_true) != len(y_pred):
        raise ValueError(f"y_true and y_pred must have the same length, got {len(y_true)} and {len(y_pred)}.")

    label_set = set(label_order)
    bad_true = sorted({lbl for lbl in y_true if lbl not in label_set})
    bad_pred = sorted({lbl for lbl in y_pred if lbl not in label_set})

    if bad_true or bad_pred:
        parts = ["Labels found that are not in the ontology (cl_names)."]
        if bad_true:
            parts.append(f"  Ground-truth labels not in ontology ({len(bad_true)} unique): {bad_true[:10]}")
        if bad_pred:
            parts.append(f"  Predicted labels not in ontology ({len(bad_pred)} unique): {bad_pred[:10]}")
        parts.append("For Azimuth outputs, curate the crosswalk so every Azimuth label maps to a valid CL term ID.")
        raise ValueError("\n".join(parts))

    dense = _sklearn_confusion_matrix(y_true, y_pred, labels=label_order)
    return scipy.sparse.csr_matrix(dense)


def save_confusion_matrix(
    cm: scipy.sparse.csr_matrix,
    meta: t.Dict[str, t.Any],
    output_dir: t.Union[str, Path],
) -> None:
    """
    Save a sparse confusion matrix and its metadata to *output_dir*.

    Writes two files:

    - ``confusion_matrix.npz`` — the sparse matrix in SciPy's NPZ format.
    - ``confusion_matrix_meta.json`` — JSON metadata (model_name, test_sample, label_order, …).

    :param cm: Sparse confusion matrix to save.
    :param meta: Metadata dict (must contain at minimum ``model_name``, ``label_order``,
        ``matrix_shape``).
    :param output_dir: Directory to write files into (created if absent).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    scipy.sparse.save_npz(str(output_dir / CM_MATRIX_FILENAME), cm.tocsr())
    with open(output_dir / CM_META_FILENAME, "w") as fh:
        json.dump(meta, fh, indent=2)


def load_confusion_matrix(
    source_dir: t.Union[str, Path],
) -> t.Tuple[scipy.sparse.csr_matrix, t.Dict[str, t.Any]]:
    """
    Load a sparse confusion matrix and its metadata from *source_dir*.

    :param source_dir: Directory containing ``confusion_matrix.npz`` and
        ``confusion_matrix_meta.json``.
    :returns: ``(cm, meta)`` tuple.
    :raises FileNotFoundError: If either required file is absent.
    """
    source_dir = Path(source_dir)
    matrix_path = source_dir / CM_MATRIX_FILENAME
    meta_path = source_dir / CM_META_FILENAME

    if not matrix_path.exists():
        raise FileNotFoundError(f"Confusion matrix not found: {matrix_path}. Run the confusion-matrix step first.")
    if not meta_path.exists():
        raise FileNotFoundError(f"Confusion matrix metadata not found: {meta_path}.")

    cm = scipy.sparse.load_npz(str(matrix_path)).tocsr()
    with open(meta_path) as fh:
        meta = json.load(fh)
    return cm, meta


def aggregate_confusion_matrices(
    cms: t.List[scipy.sparse.csr_matrix],
) -> scipy.sparse.csr_matrix:
    """
    Element-wise sum a list of confusion matrices.

    All matrices must have the same shape.  Use the ``label_order`` stored in each
    matrix's metadata to verify alignment before calling this function.

    :param cms: Non-empty list of sparse confusion matrices with identical shapes.
    :raises ValueError: If *cms* is empty or matrices have mismatched shapes.
    :returns: Summed CSR sparse matrix.
    """
    if not cms:
        raise ValueError("Cannot aggregate an empty list of confusion matrices.")

    shapes = {cm.shape for cm in cms}
    if len(shapes) > 1:
        raise ValueError(f"All confusion matrices must have the same shape; found shapes: {shapes}.")

    result = cms[0].astype(np.int64)
    for cm in cms[1:]:
        result = result + cm.astype(np.int64)
    return result.tocsr()
