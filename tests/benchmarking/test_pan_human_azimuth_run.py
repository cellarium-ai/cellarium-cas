"""Unit tests for Pan-Human Azimuth runner preprocessing, native validation, and writing."""

import json

import anndata
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellarium.cas.benchmarking.comparators.pan_human_azimuth import run


def _adata(*, x, var_index, obs_index=None, var=None, counts=None):
    x = np.asarray(x, dtype=float)
    obs_index = obs_index or [f"cell{i}" for i in range(x.shape[0])]
    adata = anndata.AnnData(
        X=x, obs=pd.DataFrame(index=obs_index), var=pd.DataFrame(var or {}, index=pd.Index(var_index))
    )
    if counts is not None:
        adata.layers["counts"] = counts
    return adata


def _native(*, barcodes, fine="Fine", medium="Medium", broad="Broad"):
    return pd.DataFrame(
        {
            "azimuth_fine": fine,
            "azimuth_medium": medium,
            "azimuth_broad": broad,
            "full_hierarchical_labels": "Broad|Medium|Fine",
            "final_level_labels": "Medium",
            "final_level_confidence": 0.9,
            "full_consistent_hierarchy": True,
        },
        index=pd.Index(barcodes),
    )


@pytest.fixture
def fake_model(monkeypatch):
    """Replace only the external model boundary; preprocessing and validation stay real."""
    monkeypatch.setattr(run, "_model_genes", lambda: {"A", "B"})
    monkeypatch.setattr(run, "_api_version", lambda: "9.9.9")
    state = {"native": lambda prepared: _native(barcodes=prepared.obs_names), "calls": 0}

    def annotate(*, prepared, batch_size):
        state["calls"] += 1
        return state["native"](prepared)

    monkeypatch.setattr(run, "_annotate", annotate)
    return state


def _write_h5ad(path, **kwargs):
    _adata(**kwargs).write_h5ad(path)
    return str(path)


@pytest.mark.parametrize("to_sparse", [False, True])
def test_counts_layer_wins_and_duplicate_symbols_sum_before_normalization(to_sparse):
    counts = np.array([[1, 3, 6], [0, 0, 0]], dtype=float)
    adata = _adata(
        x=[[99, 99, 99], [99, 99, 99]],
        var_index=["id1", "id2", "id3"],
        obs_index=["bc1", "bc2"],
        var={"feature_name": ["A", "A", "B"]},
        counts=sparse.csr_matrix(counts) if to_sparse else counts,
    )

    prepared = run._prepare_dataset(adata=adata)

    assert prepared.obs_names.tolist() == ["bc1", "bc2"]
    assert prepared.var["feature_name"].tolist() == ["A", "B"]
    np.testing.assert_allclose(prepared.X.toarray(), [[np.log1p(4000), np.log1p(6000)], [0, 0]], rtol=1e-6)


def test_x_used_when_counts_layer_absent_and_dense_integer_floats_accepted():
    adata = _adata(x=[[2.0, 8.0]], var_index=["A", "B"])

    prepared = run._prepare_dataset(adata=adata)

    np.testing.assert_allclose(prepared.X.toarray(), [[np.log1p(2000), np.log1p(8000)]], rtol=1e-6)


@pytest.mark.parametrize("bad", [0.5, -1.0, np.nan, np.inf])
def test_non_raw_counts_rejected(bad):
    adata = _adata(x=[[1.0, bad]], var_index=["A", "B"])

    with pytest.raises(ValueError, match="Raw counts"):
        run._prepare_dataset(adata=adata)


def test_symbol_precedence_and_blank_symbols():
    var = {"feature_name": ["feat_a", "feat_b"], "gene_symbol": ["sym_a", "sym_b"]}
    adata = _adata(x=[[1, 1]], var_index=["id_a", "id_b"], var=var)

    assert run._gene_symbols(adata=adata).tolist() == ["feat_a", "feat_b"]
    del adata.var["feature_name"]
    assert run._gene_symbols(adata=adata).tolist() == ["sym_a", "sym_b"]
    del adata.var["gene_symbol"]
    assert run._gene_symbols(adata=adata).tolist() == ["id_a", "id_b"]

    for blank in ("", " ", None):
        bad = _adata(x=[[1, 1]], var_index=["a", "b"], var={"feature_name": ["A", blank]})
        with pytest.raises(ValueError, match="missing or empty"):
            run._gene_symbols(adata=bad)


def test_empty_and_duplicate_barcode_inputs_rejected():
    with pytest.raises(ValueError, match="no cells"):
        run._prepare_dataset(adata=_adata(x=np.zeros((0, 2)), var_index=["A", "B"]))
    with pytest.raises(ValueError, match="no genes"):
        run._prepare_dataset(adata=_adata(x=np.zeros((2, 0)), var_index=[]))
    with pytest.raises(ValueError, match="duplicate cell barcodes"):
        run._prepare_dataset(adata=_adata(x=[[1], [1]], var_index=["A"], obs_index=["x", "x"]))


def test_unassigned_in_any_level_rejects_dataset_with_bounded_unique_diagnostics():
    barcodes = pd.Index([f"bc{i}" for i in range(8)])
    native = _native(barcodes=barcodes)
    native.loc[barcodes[:7], "azimuth_fine"] = "Unassigned"
    native.loc[barcodes[:2], "azimuth_broad"] = "Unassigned"  # same cells, counted once

    with pytest.raises(ValueError) as error:
        run._validate_native(cells_meta=native, barcodes=barcodes, input_h5ad="in.h5ad")

    assert str(error.value) == (
        "Pan-Human Azimuth produced Unassigned annotations for 7 cells in 'in.h5ad' "
        "(first 5 barcodes): ['bc0', 'bc1', 'bc2', 'bc3', 'bc4']"
    )
    only_broad = _native(barcodes=barcodes, broad="Unassigned")
    with pytest.raises(ValueError, match="for 8 cells"):
        run._validate_native(cells_meta=only_broad, barcodes=barcodes, input_h5ad="in.h5ad")


def test_validate_native_reorders_and_checks_barcodes_and_columns():
    barcodes = pd.Index(["a", "b", "c"])
    native = _native(barcodes=barcodes)
    native["azimuth_fine"] = ["fa", "fb", "fc"]

    aligned = run._validate_native(cells_meta=native.iloc[::-1], barcodes=barcodes, input_h5ad="x")
    assert aligned.index.tolist() == ["a", "b", "c"]
    assert aligned["azimuth_fine"].tolist() == ["fa", "fb", "fc"]

    with pytest.raises(ValueError, match="1 missing"):
        run._validate_native(cells_meta=native.iloc[:2], barcodes=barcodes, input_h5ad="x")
    extra = pd.concat([native, _native(barcodes=["d"])])
    with pytest.raises(ValueError, match="1 unexpected"):
        run._validate_native(cells_meta=extra, barcodes=barcodes, input_h5ad="x")
    dup = pd.concat([native, native.iloc[:1]])
    with pytest.raises(ValueError, match="duplicate"):
        run._validate_native(cells_meta=dup, barcodes=barcodes, input_h5ad="x")
    with pytest.raises(ValueError, match="missing columns"):
        run._validate_native(
            cells_meta=native.drop(columns="final_level_confidence"), barcodes=barcodes, input_h5ad="x"
        )
    empty = native.copy()
    empty["azimuth_medium"] = ["m", "", "m"]
    with pytest.raises(ValueError, match="azimuth_medium"):
        run._validate_native(cells_meta=empty, barcodes=barcodes, input_h5ad="x")


def test_run_writes_native_columns_unchanged_in_input_order(tmp_path, fake_model):
    # Native final label differs from refined fine; repeated label across levels.
    def native(prepared):
        meta = _native(barcodes=prepared.obs_names, fine="T cell", medium="T cell", broad="Lymphocyte")
        meta["final_level_labels"] = "CD4 T"
        meta["final_level_confidence"] = np.linspace(0.5, 0.75, len(meta))
        return meta.iloc[::-1]

    fake_model["native"] = native
    path = _write_h5ad(tmp_path / "s.h5ad", x=[[1, 2, 3], [4, 5, 6]], var_index=["A", "B", "C"], obs_index=["z", "y"])

    stats = run.run_pan_human_azimuth(input_h5ad=(path,), output_dir_template=str(tmp_path / "out_{stem}"))

    out = tmp_path / "out_s"
    predictions = pd.read_csv(out / run.PREDICTIONS_FILENAME, index_col=0)
    assert predictions.columns.tolist() == list(run.NATIVE_COLUMNS)
    assert predictions.index.tolist() == ["z", "y"]
    assert predictions["azimuth_fine"].tolist() == ["T cell", "T cell"]
    assert predictions["final_level_labels"].tolist() == ["CD4 T", "CD4 T"]
    np.testing.assert_allclose(predictions["final_level_confidence"], [0.5, 0.75])
    written = json.loads((out / run.STATS_FILENAME).read_text())
    assert written == {k: v for k, v in stats[0].items() if k != "stem"}
    assert written["n_cells"] == 2
    assert written["n_input_genes"] == 3
    assert written["n_model_genes"] == 2
    assert written["n_matched_genes"] == 2
    assert written["overlap_ratio"] == 1.0
    assert (written["api_version"], written["model_version"]) == ("9.9.9", "v1")


def test_zero_gene_overlap_fails_before_inference(tmp_path, fake_model):
    path = _write_h5ad(tmp_path / "s.h5ad", x=[[1, 2]], var_index=["X", "Y"])

    with pytest.raises(ValueError, match="overlap"):
        run.run_pan_human_azimuth(input_h5ad=(path,), output_dir_template=str(tmp_path / "o_{stem}"))

    assert fake_model["calls"] == 0
    assert not (tmp_path / "o_s").exists()


def test_invalid_options_fail_before_any_work(tmp_path, fake_model):
    path = _write_h5ad(tmp_path / "s.h5ad", x=[[1, 2]], var_index=["A", "B"])
    template = str(tmp_path / "o_{stem}")

    for kwargs in (
        {"input_h5ad": (), "output_dir_template": template},
        {"input_h5ad": (path,), "output_dir_template": str(tmp_path / "o")},
        {"input_h5ad": (str(tmp_path / "missing.h5ad"),), "output_dir_template": template},
        {"input_h5ad": (path,), "output_dir_template": template, "batch_size": 0},
    ):
        with pytest.raises(ValueError):
            run.run_pan_human_azimuth(**kwargs)
    assert fake_model["calls"] == 0


def test_output_collision_and_nonempty_directory_fail_without_overwrite(tmp_path, fake_model):
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    first = _write_h5ad(tmp_path / "a" / "s.h5ad", x=[[1, 2]], var_index=["A", "B"])
    second = _write_h5ad(tmp_path / "b" / "s.h5ad", x=[[1, 2]], var_index=["A", "B"])

    with pytest.raises(ValueError, match="same output directory"):
        run.run_pan_human_azimuth(input_h5ad=(first, second), output_dir_template=str(tmp_path / "o_{stem}"))

    existing = tmp_path / "o_s"
    existing.mkdir()
    (existing / "keep.txt").write_text("keep")
    with pytest.raises(ValueError, match="not empty"):
        run.run_pan_human_azimuth(input_h5ad=(first,), output_dir_template=str(tmp_path / "o_{stem}"))

    assert (existing / "keep.txt").read_text() == "keep"
    assert [p.name for p in existing.iterdir()] == ["keep.txt"]
    assert fake_model["calls"] == 0


def test_later_unassigned_input_leaves_earlier_outputs_and_writes_nothing_for_itself(tmp_path, fake_model):
    good = _write_h5ad(tmp_path / "good.h5ad", x=[[1, 2]], var_index=["A", "B"], obs_index=["g"])
    bad = _write_h5ad(tmp_path / "bad.h5ad", x=[[1, 2]], var_index=["A", "B"], obs_index=["b"])

    def native(prepared):
        return _native(barcodes=prepared.obs_names, medium="Unassigned" if prepared.obs_names[0] == "b" else "ok")

    fake_model["native"] = native

    with pytest.raises(ValueError, match="Unassigned annotations for 1 cells"):
        run.run_pan_human_azimuth(input_h5ad=(good, bad), output_dir_template=str(tmp_path / "o_{stem}"))

    assert (tmp_path / "o_good" / run.PREDICTIONS_FILENAME).is_file()
    assert (tmp_path / "o_good" / run.STATS_FILENAME).is_file()
    assert not (tmp_path / "o_bad").exists()
