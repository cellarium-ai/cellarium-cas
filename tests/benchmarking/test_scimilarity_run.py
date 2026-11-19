"""Unit tests for SCimilarity runner preprocessing and ranked output."""

import json

import anndata
import numpy as np
import pandas as pd
import pytest

from cellarium.cas.benchmarking.comparators.scimilarity import run


def _adata(*, x, var_index, var=None):
    """Build minimal AnnData for runner tests."""
    return anndata.AnnData(X=np.asarray(x), var=pd.DataFrame(var or {}, index=var_index))


def test_raw_counts_prefers_counts_layer():
    adata = _adata(x=[[99.0]], var_index=["GENE1"])
    adata.layers["counts"] = np.array([[2.0]])

    np.testing.assert_array_equal(run._raw_counts(adata=adata), [[2.0]])


def test_gene_symbols_follow_documented_precedence():
    adata = _adata(
        x=[[1.0, 2.0]],
        var_index=["var_a", "var_b"],
        var={"feature_name": ["feature_a", "feature_b"], "gene_symbol": ["symbol_a", "symbol_b"]},
    )

    assert run._gene_symbols(adata=adata).tolist() == ["feature_a", "feature_b"]
    del adata.var["feature_name"]
    assert run._gene_symbols(adata=adata).tolist() == ["symbol_a", "symbol_b"]
    del adata.var["gene_symbol"]
    assert run._gene_symbols(adata=adata).tolist() == ["var_a", "var_b"]


def test_prepare_passes_duplicate_symbols_to_consolidation(monkeypatch):
    adata = _adata(
        x=[[1.0, 2.0]],
        var_index=["id_a", "id_b"],
        var={"feature_name": ["GENE1", "GENE1"]},
    )
    captured = {}

    def consolidate_duplicate_symbols(*, adata):
        captured["symbols"] = adata.var_names.tolist()
        assert "counts" in adata.layers
        return adata

    monkeypatch.setattr(
        run,
        "_load_scimilarity_utilities",
        lambda: (lambda **kwargs: kwargs["data"], consolidate_duplicate_symbols, lambda **kwargs: kwargs["data"]),
    )

    run._prepare_dataset(adata=adata, gene_order=["GENE1"], gene_overlap_threshold=1)

    assert captured["symbols"] == ["GENE1", "GENE1"]


def test_prepare_rejects_insufficient_gene_overlap(monkeypatch):
    adata = _adata(x=[[1.0]], var_index=["GENE1"])
    monkeypatch.setattr(
        run,
        "_load_scimilarity_utilities",
        lambda: (lambda **kwargs: kwargs["data"], lambda **kwargs: kwargs["adata"], lambda **kwargs: kwargs["data"]),
    )

    with pytest.raises(ValueError, match="Gene overlap of 1 is less than --gene-overlap-threshold 2"):
        run._prepare_dataset(adata=adata, gene_order=["GENE1", "GENE2"], gene_overlap_threshold=2)


def test_weighted_hits_rank_stably_and_normalize():
    ranked = run._rank_weighted_hits(
        weighted_hits=json.dumps({"first": 2.0, "second": 2.0, "third": 1.0}),
        top_k=3,
    )

    assert ranked == [("first", 0.4), ("second", 0.4), ("third", 0.2)]


def test_top_rank_matches_upstream_weighted_prediction():
    predictions = pd.Series(["second"])
    stats = pd.DataFrame({"hits_weighted": [json.dumps({"second": 3.0, "first": 1.0})]})

    output = run._predictions_dataframe(
        obs_names=["cell1"],
        predictions=predictions,
        prediction_stats=stats,
        top_k=3,
    )

    assert output.loc["cell1"].tolist() == ["second", 0.75, "first", 0.25, None, None]
