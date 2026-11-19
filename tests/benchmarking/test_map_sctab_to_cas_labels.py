"""
Unit tests for mapping scTab top-k predictions to CAS inferred labels.
"""

import json

import anndata
import numpy as np
import pandas as pd
import pytest

from cellarium.cas.benchmarking.comparators.sctab.helpers.map_sctab_to_cas_labels import map_sctab_to_cas_labels

from ._fixtures import MOCK_ONTOLOGY_RESOURCE

BARCODES = ["cell1", "cell2", "cell3"]
TOP_K = 3

# Crosswalk: two distinct scTab labels, both mapping into the mock ontology.
CROSSWALK_ROWS = [
    {"label": "alpha", "cl_id": "CL:0000002", "cl_label": "child_a"},
    {"label": "beta", "cl_id": "CL:0000004", "cl_label": "grandchild"},
]


def _write_ontology_resource(tmp_path):
    path = tmp_path / "ontology_resource.json"
    with open(path, "w") as f:
        json.dump(MOCK_ONTOLOGY_RESOURCE, f)
    return str(path)


def _write_h5ad(tmp_path):
    path = tmp_path / "cells.h5ad"
    adata = anndata.AnnData(
        X=np.zeros((len(BARCODES), 1)),
        obs=pd.DataFrame(index=BARCODES),
        var=pd.DataFrame(index=["gene1"]),
    )
    adata.write_h5ad(str(path))
    return str(path)


def _write_predictions(tmp_path, rows):
    path = tmp_path / "predictions.csv"
    df = pd.DataFrame(rows, index=pd.Index(BARCODES, name=None))
    df.to_csv(str(path))
    return str(path)


def _write_crosswalk(tmp_path, rows=CROSSWALK_ROWS):
    path = tmp_path / "crosswalk.csv"
    pd.DataFrame(rows).to_csv(str(path), index=False)
    return str(path)


def _call_map(tmp_path, predictions_rows, crosswalk_rows=CROSSWALK_ROWS, model_name="sctab_tabnet"):
    return map_sctab_to_cas_labels(
        predictions_csv_path=_write_predictions(tmp_path, predictions_rows),
        h5ad_path=_write_h5ad(tmp_path),
        output_dir=str(tmp_path / "out"),
        crosswalk_csv_path=_write_crosswalk(tmp_path, crosswalk_rows),
        crosswalk_label_col="label",
        crosswalk_cl_id_col="cl_id",
        ontology_resource_path=_write_ontology_resource(tmp_path),
        model_name=model_name,
        crosswalk_cl_label_col="cl_label",
        top_k=TOP_K,
    )


def _rank_row(rank1_label, rank1_score, rank2_label, rank2_score, rank3_label, rank3_score):
    row = {f"predicted_label_{k}": label for k, label in [(1, rank1_label), (2, rank2_label), (3, rank3_label)]}
    row.update({f"predicted_score_{k}": score for k, score in [(1, rank1_score), (2, rank2_score), (3, rank3_score)]})
    return row


def test_unmapped_label_raises(tmp_path):
    rows = [
        _rank_row("not_in_crosswalk", 0.9, "alpha", 0.05, "beta", 0.05),
        _rank_row("alpha", 0.8, "beta", 0.1, "alpha", 0.1),
        _rank_row("beta", 0.7, "alpha", 0.2, "beta", 0.1),
    ]
    with pytest.raises(ValueError, match="missing from the crosswalk"):
        _call_map(tmp_path, rows)


def test_non_ontology_cl_id_raises(tmp_path):
    rows = [
        _rank_row("alpha", 0.9, "beta", 0.05, "alpha", 0.05),
        _rank_row("alpha", 0.8, "beta", 0.1, "alpha", 0.1),
        _rank_row("beta", 0.7, "alpha", 0.2, "beta", 0.1),
    ]
    crosswalk_rows = CROSSWALK_ROWS + [{"label": "gamma", "cl_id": "CL:9999999", "cl_label": "nope"}]
    bad_rows = [_rank_row("gamma", 0.6, "alpha", 0.3, "beta", 0.1)] + rows[:2]
    with pytest.raises(ValueError, match="absent from the ontology resource"):
        _call_map(tmp_path, bad_rows, crosswalk_rows=crosswalk_rows)


def test_duplicate_cl_ids_sum_and_rerank(tmp_path):
    # alpha and beta both map to the same CL id via the crosswalk rows below; their
    # per-rank scores must sum into one output rank, ordered by merged score.
    crosswalk_rows = [
        {"label": "alpha", "cl_id": "CL:0000002", "cl_label": "child_a"},
        {"label": "beta", "cl_id": "CL:0000002", "cl_label": "child_a"},
    ]
    rows = [
        _rank_row("alpha", 0.6, "beta", 0.3, "beta", 0.1),
        _rank_row("alpha", 0.6, "beta", 0.3, "beta", 0.1),
        _rank_row("alpha", 0.6, "beta", 0.3, "beta", 0.1),
    ]
    result = _call_map(tmp_path, rows, crosswalk_rows=crosswalk_rows)
    labels_df = pd.read_csv(result["inferred_labels_path"], index_col=0)
    assert labels_df.at["cell1", "cas_cell_type_name_1"] == "CL:0000002"
    assert labels_df.at["cell1", "cas_cell_type_score_1"] == pytest.approx(1.0)
    # Only one distinct CL id survives; rank 2 is empty.
    assert pd.isna(labels_df.at["cell1", "cas_cell_type_name_2"])


def test_metadata_model_name_and_rank_columns(tmp_path):
    rows = [
        _rank_row("alpha", 0.9, "beta", 0.05, "beta", 0.05),
        _rank_row("alpha", 0.8, "beta", 0.1, "alpha", 0.1),
        _rank_row("beta", 0.7, "alpha", 0.2, "beta", 0.1),
    ]
    result = _call_map(tmp_path, rows)
    with open(result["metadata_path"]) as f:
        metadata = json.load(f)
    assert metadata["model_name"] == "sctab_tabnet"
    assert metadata["n_cells"] == 3

    labels_df = pd.read_csv(result["inferred_labels_path"], index_col=0)
    assert list(labels_df.columns) == [
        "cas_cell_type_label_1",
        "cas_cell_type_label_2",
        "cas_cell_type_label_3",
        "cas_cell_type_name_1",
        "cas_cell_type_name_2",
        "cas_cell_type_name_3",
        "cas_cell_type_score_1",
        "cas_cell_type_score_2",
        "cas_cell_type_score_3",
    ]
    # Distinct top-1 CL id and score carried through unchanged (no merge on rank 1).
    assert labels_df.at["cell1", "cas_cell_type_name_1"] == "CL:0000002"
    assert labels_df.at["cell1", "cas_cell_type_score_1"] == pytest.approx(0.9)
    # cl_label from the crosswalk lands in the label column.
    assert labels_df.at["cell1", "cas_cell_type_label_1"] == "child_a"
