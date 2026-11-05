"""Unit tests for mapping Pan-Human Azimuth predictions to CAS labels."""

import json

import anndata
import numpy as np
import pandas as pd
import pytest

from cellarium.cas.benchmarking.comparators.pan_human_azimuth.helpers.pan_human_azimuth_to_cas_annotation import (
    pan_human_azimuth_to_cas_annotation,
)
from tests.benchmarking._fixtures import MOCK_ONTOLOGY_RESOURCE


def _convert(tmp_path, *, predictions, crosswalk):
    h5ad = tmp_path / "cells.h5ad"
    anndata.AnnData(X=np.ones((2, 1)), obs=pd.DataFrame(index=["cell2", "cell1"])).write_h5ad(h5ad)
    pd.DataFrame(predictions).to_csv(tmp_path / "p.csv")
    pd.DataFrame(crosswalk).to_csv(tmp_path / "x.csv", index=False)
    (tmp_path / "o.json").write_text(json.dumps(MOCK_ONTOLOGY_RESOURCE))
    return pan_human_azimuth_to_cas_annotation(
        predictions_csv_path=str(tmp_path / "p.csv"),
        h5ad_path=str(h5ad),
        crosswalk_csv_path=str(tmp_path / "x.csv"),
        crosswalk_label_col="pan_human_label",
        crosswalk_cl_id_col="cl_id",
        crosswalk_cl_label_col="cl_label",
        ontology_resource_path=str(tmp_path / "o.json"),
        output_dir=str(tmp_path / "out"),
    )


CROSSWALK = [
    {"pan_human_label": "fa", "cl_id": "CL:0000002", "cl_label": "child_a"},
    {"pan_human_label": "fb", "cl_id": "CL:0000003", "cl_label": "child_b"},
    {"pan_human_label": "mid", "cl_id": "CL:0000002", "cl_label": "child_a"},
    {"pan_human_label": "top", "cl_id": "CL:0000003", "cl_label": "child_b"},
]


def test_levels_become_ranks_fine_first_with_blank_scores_in_h5ad_order(tmp_path):
    predictions = {
        "azimuth_fine": {"cell1": "fa", "cell2": "fb"},
        "azimuth_medium": {"cell1": "mid", "cell2": "mid"},
        "azimuth_broad": {"cell1": "top", "cell2": "top"},
        "final_level_labels": {"cell1": "ignored", "cell2": "ignored"},
    }

    result = _convert(tmp_path, predictions=predictions, crosswalk=CROSSWALK)

    labels = pd.read_csv(result["inferred_labels_path"], index_col=0)
    assert labels.index.tolist() == ["cell2", "cell1"]
    assert labels.columns.tolist() == [
        f"cas_cell_type_{kind}_{rank}" for kind in ("label", "name", "score") for rank in (1, 2, 3)
    ]
    # Same CL term at two ranks is kept per level, not merged.
    assert labels.loc["cell1", ["cas_cell_type_name_1", "cas_cell_type_name_2", "cas_cell_type_name_3"]].tolist() == [
        "CL:0000002",
        "CL:0000002",
        "CL:0000003",
    ]
    assert labels.loc["cell2", "cas_cell_type_name_1"] == "CL:0000003"
    assert labels[[f"cas_cell_type_score_{k}" for k in (1, 2, 3)]].isna().all().all()
    assert json.loads((tmp_path / "out" / "metadata.json").read_text())["model_name"] == "pan_human_azimuth_v1"
    assert (tmp_path / "out" / "ontology_response.json").is_file()


def test_rejects_unmapped_labels_in_any_level_and_missing_columns(tmp_path):
    predictions = {
        "azimuth_fine": {"cell1": "fa", "cell2": "fb"},
        "azimuth_medium": {"cell1": "mid", "cell2": "mid"},
        "azimuth_broad": {"cell1": "top", "cell2": "unknown"},
    }
    with pytest.raises(ValueError, match=r"Pan-Human Azimuth class labels missing from the crosswalk: \['unknown'\]"):
        _convert(tmp_path, predictions=predictions, crosswalk=CROSSWALK)

    with pytest.raises(ValueError, match="missing columns"):
        _convert(tmp_path, predictions={"azimuth_fine": {"cell1": "fa", "cell2": "fb"}}, crosswalk=CROSSWALK)
