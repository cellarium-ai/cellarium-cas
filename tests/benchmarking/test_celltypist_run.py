"""Unit tests for CellTypist runner preprocessing, ranking, and option validation."""

import anndata
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellarium.cas.benchmarking.comparators.celltypist import run


def _adata(*, x, var_index, obs_index=None, counts=None):
    x = np.asarray(x, dtype=float)
    obs_index = obs_index or [f"cell{i}" for i in range(x.shape[0])]
    adata = anndata.AnnData(X=x, obs=pd.DataFrame(index=obs_index), var=pd.DataFrame(index=pd.Index(var_index)))
    if counts is not None:
        adata.layers["counts"] = counts
    return adata


@pytest.mark.parametrize("to_sparse", [False, True])
def test_prepare_uses_counts_layer_sums_duplicates_and_normalizes(to_sparse):
    counts = np.array([[1, 1, 2, 0], [0, 0, 0, 0], [3, 0, 0, 1]], dtype=float)
    layer = sparse.csr_matrix(counts) if to_sparse else counts
    adata = _adata(x=np.full((3, 4), 0.5), var_index=["A", "A", "B", "C"], obs_index=["z", "y", "x"], counts=layer)

    prepared = run._prepare_dataset(adata=adata)

    assert prepared.obs_names.tolist() == ["z", "y", "x"]
    assert prepared.var_names.tolist() == ["A", "B", "C"]
    dense = np.expm1(prepared.X.toarray() if sparse.issparse(prepared.X) else prepared.X)
    np.testing.assert_allclose(dense[0], [2 / 4 * 1e4, 2 / 4 * 1e4, 0], rtol=1e-5)
    np.testing.assert_allclose(dense.sum(axis=1)[[0, 2]], 1e4, rtol=1e-5)
    assert dense[1].sum() == 0  # zero-total cell retained


@pytest.mark.parametrize(
    "x, obs, message",
    [
        (np.array([[1.0, np.nan]]), ["c"], "finite"),
        (np.array([[1.0, -1.0]]), ["c"], "negative"),
        (np.array([[1.0, 0.5]]), ["c"], "integer-valued"),
        (np.array([[1.0, 1.0], [1.0, 1.0]]), ["c", "c"], "duplicate cell barcodes"),
    ],
)
def test_prepare_rejects_invalid_input(x, obs, message):
    with pytest.raises(ValueError, match=message):
        run._prepare_dataset(adata=_adata(x=x, var_index=["A", "B"], obs_index=obs))


def _matrices(decision, probability, *, index, classes):
    return (
        pd.DataFrame({"predicted_labels": [classes[i] for i in np.argmax(decision, axis=1)]}, index=index),
        pd.DataFrame(decision, index=index, columns=classes),
        pd.DataFrame(probability, index=index, columns=classes),
    )


def test_rank_native_orders_by_decision_not_saturated_probability_and_pads():
    index = ["c2", "c1"]
    classes = ["A", "B", "C"]
    # Row c2: probabilities saturate to 1.0 for all, decision order is B > C > A.
    # Row c1: exact decision tie between A and C; model-column order must win.
    decision = np.array([[-20.0, 30.0, 25.0], [4.0, -1.0, 4.0]])
    probability = np.array([[1.0, 1.0, 1.0], [0.98, 0.27, 0.98]])
    labels, dec, prob = _matrices(decision, probability, index=index, classes=classes)

    aligned, ranked = run._rank_native(
        predicted_labels=labels.loc[["c2", "c1"]],
        decision=dec.loc[["c1", "c2"]],  # shuffled upstream order
        probability=prob.loc[["c1", "c2"]],
        barcodes=pd.Index(["c2", "c1"]),
        top_k=4,
    )

    assert ranked.index.tolist() == ["c2", "c1"]
    assert aligned.index.tolist() == ["c2", "c1"]
    assert ranked.loc["c2", ["predicted_label_1", "predicted_label_2", "predicted_label_3"]].tolist() == ["B", "C", "A"]
    assert ranked.loc["c1", ["predicted_label_1", "predicted_label_2", "predicted_label_3"]].tolist() == ["A", "C", "B"]
    assert ranked.at["c1", "predicted_score_3"] == pytest.approx(0.27)
    assert pd.isna(ranked.at["c1", "predicted_label_4"]) and pd.isna(ranked.at["c1", "predicted_score_4"])


@pytest.mark.parametrize(
    "mutate, message",
    [
        (lambda p, d, l: (p, d, l.iloc[:1]), "one-to-one"),
        (lambda p, d, l: (p.assign(A=np.nan), d, l), "non-finite"),
        (lambda p, d, l: (p.assign(A=1.5), d, l), r"\[0, 1\]"),
        (lambda p, d, l: (p, d, l.assign(predicted_labels="B")), "disagrees"),
    ],
)
def test_rank_native_rejects_invalid_results(mutate, message):
    index, classes = ["c1", "c2"], ["A", "B"]
    labels, dec, prob = _matrices(
        np.array([[2.0, 1.0], [2.0, 1.0]]), np.array([[0.9, 0.7], [0.9, 0.7]]), index=index, classes=classes
    )
    prob, dec, labels = mutate(prob, dec, labels)
    with pytest.raises(ValueError, match=message):
        run._rank_native(predicted_labels=labels, decision=dec, probability=prob, barcodes=pd.Index(index), top_k=2)


def test_run_rejects_bad_options_before_loading_model(tmp_path):
    model = tmp_path / "m.pkl"
    model.write_bytes(b"x")
    h5ad = tmp_path / "s.h5ad"
    _adata(x=[[1.0]], var_index=["A"]).write_h5ad(h5ad)
    occupied = tmp_path / "out_s"
    occupied.mkdir()
    (occupied / "keep").write_text("1")

    with pytest.raises(ValueError, match="top_k"):
        run.run_celltypist(
            model_path=str(model), input_h5ad=(str(h5ad),), output_dir_template=str(tmp_path / "{stem}"), top_k=0
        )
    with pytest.raises(ValueError, match="Model file"):
        run.run_celltypist(
            model_path=str(tmp_path / "no.pkl"), input_h5ad=(str(h5ad),), output_dir_template=str(tmp_path / "{stem}")
        )
    with pytest.raises(ValueError, match=r"\{stem\}"):
        run.run_celltypist(model_path=str(model), input_h5ad=(str(h5ad),), output_dir_template=str(tmp_path / "out"))
    with pytest.raises(ValueError, match="not empty"):
        run.run_celltypist(
            model_path=str(model), input_h5ad=(str(h5ad),), output_dir_template=str(tmp_path / "out_{stem}")
        )
