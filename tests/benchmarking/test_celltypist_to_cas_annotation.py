"""Unit tests for converting native CellTypist predictions to CAS annotations."""

import json

import anndata
import numpy as np
import pandas as pd
import pytest

from cellarium.cas.benchmarking.comparators.celltypist.helpers.celltypist_to_cas_annotation import (
    celltypist_to_cas_annotation,
)
from tests.benchmarking._fixtures import MOCK_ONTOLOGY_RESOURCE

CROSSWALK = [
    {"celltypist_label": "A", "cl_id": "CL:0000002", "cl_label": "child_a"},
    {"celltypist_label": "B", "cl_id": "CL:0000003", "cl_label": "child_b"},
    {"celltypist_label": "C", "cl_id": "CL:0000003", "cl_label": "child_b"},
]


def _setup(tmp_path, *, crosswalk=CROSSWALK, stats_update=None, predictions=None):
    h5ad = tmp_path / "cells.h5ad"
    anndata.AnnData(X=np.ones((2, 1)), obs=pd.DataFrame(index=["cell2", "cell1"])).write_h5ad(h5ad)
    native = tmp_path / "native"
    native.mkdir()
    if predictions is None:
        predictions = pd.DataFrame(  # barcode order reversed relative to the H5AD
            {
                "predicted_label_1": ["A", "B"],
                "predicted_score_1": [0.8, 0.9],
                "predicted_label_2": ["B", "C"],
                "predicted_score_2": [0.6, 0.5],
                "predicted_label_3": ["C", None],
                "predicted_score_3": [0.5, None],
            },
            index=["cell1", "cell2"],
        )
    predictions.to_csv(native / "celltypist_predictions.csv")
    stats = {
        "input_h5ad": str(h5ad.resolve()),
        "model_name": "celltypist_m_abc_raw",
        "model_sha256": "abc",
        "api_version": "1.7.1",
        "mode": "best match",
        "majority_voting": False,
        "top_k": 3,
    }
    stats.update(stats_update or {})
    (native / "celltypist_stats.json").write_text(json.dumps(stats))
    crosswalk_path = tmp_path / "crosswalk.csv"
    pd.DataFrame(crosswalk).to_csv(crosswalk_path, index=False)
    ontology = tmp_path / "ontology_resource.json"
    ontology.write_text(json.dumps(MOCK_ONTOLOGY_RESOURCE))
    return dict(
        native_output_dir=str(native),
        h5ad_path=str(h5ad),
        crosswalk_csv_path=str(crosswalk_path),
        crosswalk_label_col="celltypist_label",
        crosswalk_cl_id_col="cl_id",
        crosswalk_cl_label_col="cl_label",
        ontology_resource_path=str(ontology),
        output_dir=str(tmp_path / "out"),
        top_k=3,
    )


def test_duplicate_cl_ids_keep_max_score_and_rank_one(tmp_path):
    kwargs = _setup(tmp_path)

    result = celltypist_to_cas_annotation(**kwargs)

    labels = pd.read_csv(result["inferred_labels_path"], index_col=0)
    assert labels.index.tolist() == ["cell2", "cell1"]
    # cell1: A=0.8 -> CL:2; B=0.6 and C=0.5 -> CL:3 merged by max (0.6), not summed (1.1).
    assert labels.loc["cell1", ["cas_cell_type_name_1", "cas_cell_type_name_2"]].tolist() == [
        "CL:0000002",
        "CL:0000003",
    ]
    assert labels.loc["cell1", ["cas_cell_type_score_1", "cas_cell_type_score_2"]].tolist() == pytest.approx([0.8, 0.6])
    assert pd.isna(labels.at["cell1", "cas_cell_type_name_3"])
    metadata = json.loads(open(result["metadata_path"]).read())
    assert metadata["model_name"] == "celltypist_m_abc_raw"
    assert (metadata["model_sha256"], metadata["api_version"], metadata["duplicate_score_aggregation"]) == (
        "abc",
        "1.7.1",
        "max",
    )
    assert set(result) == {
        "output_dir",
        "inferred_labels_path",
        "metadata_path",
        "ontology_response_path",
        "ontology_resource_path",
    }


@pytest.mark.parametrize(
    "setup_kwargs, message",
    [
        ({"crosswalk": CROSSWALK[:2]}, "missing from the crosswalk"),
        (
            {"crosswalk": CROSSWALK + [{"celltypist_label": "A", "cl_id": "CL:0000003", "cl_label": "child_b"}]},
            "multiple targets",
        ),
        ({"crosswalk": [CROSSWALK[0] | {"cl_id": "CL:9999999"}] + CROSSWALK[1:]}, "absent from the ontology"),
        ({"crosswalk": [CROSSWALK[0] | {"cl_id": ""}] + CROSSWALK[1:]}, "empty values"),
        ({"stats_update": {"input_h5ad": "/elsewhere/other.h5ad"}}, "not"),
        ({"stats_update": {"model_name": "celltypist_m_abc_voted", "majority_voting": True}}, "raw per-cell"),
        ({"stats_update": {"top_k": 2}}, "top_k must be between"),
    ],
)
def test_conversion_fails_visibly(tmp_path, setup_kwargs, message):
    kwargs = _setup(tmp_path, **setup_kwargs)
    with pytest.raises(ValueError, match=message):
        celltypist_to_cas_annotation(**kwargs)


def _full_predictions(**overrides):
    frame = pd.DataFrame(
        {
            "predicted_label_1": ["A", "B"],
            "predicted_score_1": [0.8, 0.9],
            "predicted_label_2": ["B", "C"],
            "predicted_score_2": [0.6, 0.5],
            "predicted_label_3": ["C", None],
            "predicted_score_3": [0.5, None],
        },
        index=["cell1", "cell2"],
    )
    return frame.assign(**overrides) if overrides else frame


@pytest.mark.parametrize(
    "predictions, message",
    [
        (_full_predictions().drop(columns="predicted_score_3"), "lack columns"),
        (_full_predictions().iloc[:1], "exactly cover"),
        (_full_predictions(predicted_score_1=[0.8, 1.7]), r"within \[0, 1\]"),
    ],
)
def test_conversion_rejects_bad_native_predictions(tmp_path, predictions, message):
    with pytest.raises(ValueError, match=message):
        celltypist_to_cas_annotation(**_setup(tmp_path, predictions=predictions))


def test_conversion_refuses_nonempty_output_dir(tmp_path):
    kwargs = _setup(tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    (out / "keep.txt").write_text("x")
    with pytest.raises(ValueError, match="not empty"):
        celltypist_to_cas_annotation(**kwargs)
    assert (out / "keep.txt").read_text() == "x"
