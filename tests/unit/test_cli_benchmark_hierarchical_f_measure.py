import json
import typing as t

import anndata
import pandas as pd
import pytest

from cellarium.cas.cli.benchmark import benchmark_group, run_hierarchical_f_measure_benchmark

MOCK_ONTOLOGY_RESOURCE: t.Dict[str, t.Any] = {
    "cl_names": ["CL:0000000", "CL:0000002"],
    "cell_ontology_term_id_to_cell_type": {
        "CL:0000000": "root",
        "CL:0000002": "child_a",
    },
    "children_dictionary": {
        "CL:0000000": ["CL:0000002"],
    },
}


def write_annotate_dir(
    path,
    input_path="input.h5ad",
    prediction: t.Optional[str] = "CL:0000002",
    prediction_2: t.Optional[str] = None,
    model_name="test-model",
    include_resource=True,
    include_prediction_2=False,
):
    path.mkdir(parents=True)
    labels_data = {"cas_cell_type_name_1": [prediction]}
    if include_prediction_2:
        labels_data["cas_cell_type_name_2"] = [prediction_2]
    labels_df = pd.DataFrame(labels_data, index=[path.name])
    labels_df.to_csv(path / "inferred_labels.csv")
    (path / "metadata.json").write_text(json.dumps({"input_path": input_path, "model_name": model_name}))
    if include_resource:
        (path / "ontology_resource.json").write_text(json.dumps(MOCK_ONTOLOGY_RESOURCE))


def test_run_hierarchical_f_measure_benchmark_writes_sample_and_total_outputs(tmp_path, monkeypatch):
    parent = tmp_path / "annotate_outputs"
    write_annotate_dir(parent / "sample_a", input_path="sample_a.h5ad")
    write_annotate_dir(parent / "sample_b", input_path="sample_b.h5ad", prediction=None)
    out_dir = tmp_path / "benchmark"

    def read_h5ad(_input_path):
        return anndata.AnnData(obs=pd.DataFrame({"cell_type_ontology_term_id": ["CL:0000002"]}))

    monkeypatch.setattr(anndata, "read_h5ad", read_h5ad)

    result = run_hierarchical_f_measure_benchmark(
        annotate_dirs=parent,
        gt_cl_column_name="cell_type_ontology_term_id",
        output_dir=out_dir,
        save_class_level=True,
    )

    summary_path = out_dir / "hierarchical_f_measure_summary.csv"
    class_path = out_dir / "hierarchical_f_measure_class_level.csv"
    summary = pd.read_csv(summary_path)
    class_level = pd.read_csv(class_path)

    assert result["n_summary_rows"] == 3
    assert result["n_class_level_rows"] == 3
    assert list(summary["row_type"]) == ["sample", "sample", "total"]
    assert "top_1_macro_weighted_hierarchical_f1" in summary.columns
    assert class_path.exists()

    total = summary[summary["row_type"] == "total"].iloc[0]
    assert total["annotate_dir"] == "__total__"
    assert total["input_path"] == "__all__"
    assert total["top_1_micro_hierarchical_precision"] == pytest.approx(1.0)
    assert total["top_1_micro_hierarchical_recall"] == pytest.approx(1 / 2)
    assert total["top_1_micro_hierarchical_f1"] == pytest.approx(2 / 3)

    total_class = class_level[class_level["row_type"] == "total"].iloc[0]
    assert total_class["ground_truth_class"] == "CL:0000002"
    assert total_class["tp"] == pytest.approx(2.0)
    assert total_class["fp"] == pytest.approx(0.0)
    assert total_class["fn"] == pytest.approx(2.0)


def test_run_hierarchical_f_measure_benchmark_infers_available_top_k(tmp_path, monkeypatch):
    parent = tmp_path / "annotate_outputs"
    write_annotate_dir(
        parent / "sample_a",
        input_path="sample_a.h5ad",
        prediction="CL:0000000",
        prediction_2="CL:0000002",
        include_prediction_2=True,
    )
    out_dir = tmp_path / "benchmark"

    def read_h5ad(_input_path):
        return anndata.AnnData(obs=pd.DataFrame({"cell_type_ontology_term_id": ["CL:0000002"]}))

    monkeypatch.setattr(anndata, "read_h5ad", read_h5ad)

    run_hierarchical_f_measure_benchmark(
        annotate_dirs=parent,
        gt_cl_column_name="cell_type_ontology_term_id",
        output_dir=out_dir,
    )

    summary = pd.read_csv(out_dir / "hierarchical_f_measure_summary.csv")
    sample = summary[summary["row_type"] == "sample"].iloc[0]

    assert "top_2_micro_hierarchical_f1" in summary.columns
    assert sample["top_1_micro_hierarchical_f1"] == pytest.approx(2 / 3)
    assert sample["top_2_micro_hierarchical_f1"] == pytest.approx(1.0)


def test_hierarchical_f_measure_command_does_not_expose_top_k_option():
    command = benchmark_group.commands["hierarchical-f-measure"]

    assert "top_k" not in {param.name for param in command.params}


def test_run_hierarchical_f_measure_benchmark_requires_ontology_resource(tmp_path):
    parent = tmp_path / "annotate_outputs"
    annotate_dir = parent / "model_a"
    write_annotate_dir(annotate_dir, include_resource=False)

    with pytest.raises((ValueError, FileNotFoundError), match="ontology_resource.json"):
        run_hierarchical_f_measure_benchmark(
            annotate_dirs=parent,
            gt_cl_column_name="cell_type_ontology_term_id",
            output_dir=tmp_path / "benchmark",
        )
