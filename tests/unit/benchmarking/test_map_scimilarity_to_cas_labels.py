"""Unit tests for mapping SCimilarity predictions to CAS labels."""

import json

import anndata
import numpy as np
import pandas as pd
import pytest

from cellarium.cas.benchmarking.comparators.scimilarity.helpers.map_scimilarity_to_cas_labels import (
    map_scimilarity_to_cas_labels,
)
from tests.unit.benchmarking._fixtures import MOCK_ONTOLOGY_RESOURCE


def _write_inputs(tmp_path, *, predictions, crosswalk):
    """Write minimal prediction, crosswalk, H5AD, and ontology inputs."""
    h5ad_path = tmp_path / "cells.h5ad"
    anndata.AnnData(X=np.ones((2, 1)), obs=pd.DataFrame(index=["cell2", "cell1"])).write_h5ad(h5ad_path)
    predictions_path = tmp_path / "predictions.csv"
    pd.DataFrame(predictions).to_csv(predictions_path)
    crosswalk_path = tmp_path / "crosswalk.csv"
    pd.DataFrame(crosswalk).to_csv(crosswalk_path, index=False)
    ontology_path = tmp_path / "ontology_resource.json"
    ontology_path.write_text(json.dumps(MOCK_ONTOLOGY_RESOURCE))
    return h5ad_path, predictions_path, crosswalk_path, ontology_path


def _map(tmp_path, *, predictions, crosswalk):
    """Run the SCimilarity mapper over temporary inputs."""
    h5ad_path, predictions_path, crosswalk_path, ontology_path = _write_inputs(
        tmp_path,
        predictions=predictions,
        crosswalk=crosswalk,
    )
    return map_scimilarity_to_cas_labels(
        predictions_csv_path=str(predictions_path),
        h5ad_path=str(h5ad_path),
        output_dir=str(tmp_path / "out"),
        crosswalk_csv_path=str(crosswalk_path),
        crosswalk_label_col="scimilarity_label",
        crosswalk_cl_id_col="cl_id",
        crosswalk_cl_label_col="cl_label",
        ontology_resource_path=str(ontology_path),
        top_k=2,
    )


def test_mapper_reorders_barcodes_merges_cl_scores_and_sets_model_name(tmp_path):
    predictions = {
        "predicted_label_1": {"cell1": "alpha", "cell2": "alpha"},
        "predicted_score_1": {"cell1": 0.7, "cell2": 0.2},
        "predicted_label_2": {"cell1": "gamma", "cell2": "beta"},
        "predicted_score_2": {"cell1": 0.3, "cell2": 0.4},
    }
    crosswalk = [
        {"scimilarity_label": "alpha", "cl_id": "CL:0000002", "cl_label": "child_a"},
        {"scimilarity_label": "beta", "cl_id": "CL:0000002", "cl_label": "child_a"},
        {"scimilarity_label": "gamma", "cl_id": "CL:0000003", "cl_label": "child_b"},
    ]

    result = _map(tmp_path, predictions=predictions, crosswalk=crosswalk)

    labels = pd.read_csv(result["inferred_labels_path"], index_col=0)
    assert labels.index.tolist() == ["cell2", "cell1"]
    assert labels.columns.tolist() == [
        "cas_cell_type_label_1",
        "cas_cell_type_label_2",
        "cas_cell_type_name_1",
        "cas_cell_type_name_2",
        "cas_cell_type_score_1",
        "cas_cell_type_score_2",
    ]
    assert labels.at["cell2", "cas_cell_type_label_1"] == "child_a"
    assert pd.isna(labels.at["cell2", "cas_cell_type_label_2"])
    assert labels.at["cell2", "cas_cell_type_name_1"] == "CL:0000002"
    assert pd.isna(labels.at["cell2", "cas_cell_type_name_2"])
    assert labels.at["cell2", "cas_cell_type_score_1"] == pytest.approx(0.6)
    assert pd.isna(labels.at["cell2", "cas_cell_type_score_2"])
    assert json.loads((tmp_path / "out" / "metadata.json").read_text())["model_name"] == "scimilarity_model_v1.1"


def test_mapper_rejects_missing_crosswalk_labels(tmp_path):
    predictions = {
        "predicted_label_1": {"cell1": "alpha", "cell2": "unknown"},
        "predicted_score_1": {"cell1": 1.0, "cell2": 1.0},
        "predicted_label_2": {"cell1": None, "cell2": None},
        "predicted_score_2": {"cell1": None, "cell2": None},
    }
    crosswalk = [{"scimilarity_label": "alpha", "cl_id": "CL:0000002", "cl_label": "child_a"}]

    with pytest.raises(ValueError, match=r"SCimilarity class labels missing from the crosswalk: \['unknown'\]"):
        _map(tmp_path, predictions=predictions, crosswalk=crosswalk)
