"""Unit tests for mapping Pan-Human Azimuth predictions to CAS labels via the vendor CL map."""

import json

import anndata
import numpy as np
import pandas as pd
import pytest

from cellarium.cas.benchmarking.comparators.pan_human_azimuth.helpers.pan_human_azimuth_to_cas_annotation import (
    pan_human_azimuth_to_cas_annotation,
)
from tests.benchmarking._fixtures import MOCK_ONTOLOGY_RESOURCE

VENDOR = [
    {"Annotation_Label": "fa", "CL_ID": "CL:0000002", "CL_Label": "vendor_a"},
    {"Annotation_Label": "fb", "CL_ID": "CL:0000003", "CL_Label": "vendor_b"},
    {"Annotation_Label": "mid", "CL_ID": "CL:0000002", "CL_Label": "vendor_a"},
    {"Annotation_Label": "top", "CL_ID": "CL:0000003", "CL_Label": "vendor_b"},
]
PREDICTIONS = {
    "azimuth_fine": {"cell1": "fa", "cell2": "fb"},
    "azimuth_medium": {"cell1": "mid", "cell2": "mid"},
    "azimuth_broad": {"cell1": "top", "cell2": "top"},
}


def _convert(tmp_path, *, predictions=PREDICTIONS, vendor=VENDOR, crosswalk=None):
    h5ad = tmp_path / "cells.h5ad"
    anndata.AnnData(X=np.ones((2, 1)), obs=pd.DataFrame(index=["cell2", "cell1"])).write_h5ad(h5ad)
    pd.DataFrame(predictions).to_csv(tmp_path / "p.csv")
    pd.DataFrame(vendor).to_csv(tmp_path / "vendor.csv", index=False)
    (tmp_path / "o.json").write_text(json.dumps(MOCK_ONTOLOGY_RESOURCE))
    kwargs = {}
    if crosswalk is not None:
        pd.DataFrame(crosswalk).to_csv(tmp_path / "x.csv", index=False)
        kwargs = dict(
            crosswalk_csv_path=str(tmp_path / "x.csv"),
            crosswalk_source_cl_id_col="old_cl",
            crosswalk_target_cl_id_col="new_cl",
            crosswalk_cl_label_col="new_label",
        )
    return pan_human_azimuth_to_cas_annotation(
        predictions_csv_path=str(tmp_path / "p.csv"),
        h5ad_path=str(h5ad),
        ontology_resource_path=str(tmp_path / "o.json"),
        output_dir=str(tmp_path / "out"),
        vendor_map_csv_path=str(tmp_path / "vendor.csv"),
        **kwargs,
    )


def test_without_crosswalk_vendor_cl_ids_are_used_as_ranks_fine_first(tmp_path):
    result = _convert(tmp_path)

    labels = pd.read_csv(result["inferred_labels_path"], index_col=0)
    assert labels.index.tolist() == ["cell2", "cell1"]
    assert labels.columns.tolist() == [
        f"cas_cell_type_{kind}_{rank}" for kind in ("label", "name", "score") for rank in (1, 2, 3)
    ]
    assert labels.loc["cell1", ["cas_cell_type_name_1", "cas_cell_type_name_2", "cas_cell_type_name_3"]].tolist() == [
        "CL:0000002",
        "CL:0000002",
        "CL:0000003",
    ]
    assert labels.loc["cell2", "cas_cell_type_label_1"] == "vendor_b"
    assert labels[[f"cas_cell_type_score_{k}" for k in (1, 2, 3)]].isna().all().all()
    assert json.loads((tmp_path / "out" / "metadata.json").read_text())["model_name"] == "pan_human_azimuth_v1"
    assert (tmp_path / "out" / "ontology_response.json").is_file()


def test_crosswalk_translates_listed_cl_ids_and_passes_others_through(tmp_path):
    crosswalk = [{"old_cl": "CL:0000002", "new_cl": "CL:0000003", "new_label": "renamed"}]

    result = _convert(tmp_path, crosswalk=crosswalk)

    labels = pd.read_csv(result["inferred_labels_path"], index_col=0)
    assert labels.loc["cell1", "cas_cell_type_name_1"] == "CL:0000003"
    assert labels.loc["cell1", "cas_cell_type_label_1"] == "renamed"
    assert labels.loc["cell2", "cas_cell_type_name_1"] == "CL:0000003"
    assert labels.loc["cell2", "cas_cell_type_label_1"] == "vendor_b"


def test_rejects_labels_missing_from_vendor_map(tmp_path):
    predictions = {**PREDICTIONS, "azimuth_broad": {"cell1": "top", "cell2": "unknown"}}

    with pytest.raises(ValueError, match=r"missing from the vendor CL map: \['unknown'\]"):
        _convert(tmp_path, predictions=predictions)


def test_rejects_cl_ids_absent_from_ontology(tmp_path):
    vendor = [{**row, "CL_ID": "CL:9999999"} if row["Annotation_Label"] == "fa" else row for row in VENDOR]

    with pytest.raises(ValueError, match=r"absent from the ontology resource: \['CL:9999999'\]"):
        _convert(tmp_path, vendor=vendor)


def test_rejects_missing_prediction_columns(tmp_path):
    with pytest.raises(ValueError, match="missing columns"):
        _convert(tmp_path, predictions={"azimuth_fine": PREDICTIONS["azimuth_fine"]})
