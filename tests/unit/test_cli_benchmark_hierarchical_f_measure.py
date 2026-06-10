import json
import typing as t

import anndata
import pandas as pd
import pytest

from cellarium.cas.cli.benchmark import run_hierarchical_f_measure_benchmark
from cellarium.cas.models import CellTypeOntologyAwareResults

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


def make_match(term_id: str) -> CellTypeOntologyAwareResults.Match:
    return CellTypeOntologyAwareResults.Match(
        cell_type_ontology_term_id=term_id,
        score=1.0,
        cell_type=term_id,
    )


def make_annotation(cell_id: str, matches: t.List[CellTypeOntologyAwareResults.Match]):
    return CellTypeOntologyAwareResults.OntologyAwareAnnotation(
        query_cell_id=cell_id,
        matches=matches,
        total_weight=sum(m.score for m in matches),
        total_neighbors=len(matches),
        total_neighbors_unrecognized=0,
    )


def make_response(annotations: t.List[CellTypeOntologyAwareResults.OntologyAwareAnnotation], model_name="test-model"):
    return CellTypeOntologyAwareResults(data=annotations, model_name=model_name)


def write_annotate_dir(
    path,
    input_path="input.h5ad",
    matches: t.Sequence[str] = ("CL:0000000", "CL:0000002"),
    model_name="test-model",
    include_resource=True,
):
    path.mkdir(parents=True)
    response = make_response([make_annotation(path.name, [make_match(term_id) for term_id in matches])], model_name)
    (path / "ontology_response.json").write_text(response.model_dump_json())
    (path / "metadata.json").write_text(json.dumps({"input_path": input_path, "model_name": model_name}))
    if include_resource:
        (path / "ontology_resource.json").write_text(json.dumps(MOCK_ONTOLOGY_RESOURCE))


def test_run_hierarchical_f_measure_benchmark_writes_sample_and_total_outputs(tmp_path, monkeypatch):
    parent = tmp_path / "annotate_outputs"
    write_annotate_dir(parent / "sample_a", input_path="sample_a.h5ad")
    write_annotate_dir(parent / "sample_b", input_path="sample_b.h5ad", matches=("unknown",))
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
    assert "macro_weighted_hierarchical_f1" in summary.columns
    assert class_path.exists()

    total = summary[summary["row_type"] == "total"].iloc[0]
    assert total["annotate_dir"] == "__total__"
    assert total["input_path"] == "__all__"
    assert total["micro_hierarchical_precision"] == pytest.approx(2 / 3)
    assert total["micro_hierarchical_recall"] == pytest.approx(1 / 2)
    assert total["micro_hierarchical_f1"] == pytest.approx(4 / 7)

    total_class = class_level[class_level["row_type"] == "total"].iloc[0]
    assert total_class["ground_truth_class"] == "CL:0000002"
    assert total_class["tp"] == pytest.approx(2.0)
    assert total_class["fp"] == pytest.approx(1.0)
    assert total_class["fn"] == pytest.approx(2.0)


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
