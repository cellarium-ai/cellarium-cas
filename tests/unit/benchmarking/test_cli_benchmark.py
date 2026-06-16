import json
import typing as t
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from cellarium.cas.cli._benchmark_impl import (
    _ranked_inferred_label_columns,
    run_aggregate_step,
    run_confusion_matrix_step,
    run_f_measure_step,
    run_hierarchical_f_measure_step,
)
from cellarium.cas.cli.benchmark import benchmark_group
from tests.unit.benchmarking._fixtures import MOCK_ONTOLOGY_RESOURCE


def _write_json(path: Path, payload: t.Dict[str, t.Any]) -> None:
    with open(path, "w") as f:
        json.dump(payload, f)


def _write_benchmark_fixture(tmp_path: Path, include_rank_3: bool = True) -> t.Tuple[Path, Path]:
    input_path = tmp_path / "sample.h5ad"
    cell_ids = ["cell_1", "cell_2", "cell_3"]
    adata = anndata.AnnData(
        X=np.ones((3, 1), dtype=np.float32),
        obs=pd.DataFrame(
            {"truth": ["CL:0000002", "CL:0000003", "CL:0000004"]},
            index=cell_ids,
        ),
        var=pd.DataFrame(index=["gene_1"]),
    )
    adata.write_h5ad(input_path)

    annotate_parent = tmp_path / "annotate_outputs"
    annotate_dir = annotate_parent / "run_model_a_sample_1"
    annotate_dir.mkdir(parents=True)

    _write_json(annotate_dir / "metadata.json", {"input_path": str(input_path), "model_name": "m1"})
    _write_json(annotate_dir / "ontology_resource.json", MOCK_ONTOLOGY_RESOURCE)
    _write_json(annotate_dir / "ontology_response.json", {})

    inferred_labels = pd.DataFrame(
        {
            "cas_cell_type_name_1": ["CL:0000002", "CL:0000002", "CL:0000002"],
            "cas_cell_type_name_2": ["CL:0000001", "CL:0000003", "CL:0000004"],
        },
        index=cell_ids,
    )
    if include_rank_3:
        inferred_labels["cas_cell_type_name_3"] = ["CL:0000001", "CL:0000001", "CL:0000001"]
    inferred_labels.to_csv(annotate_dir / "inferred_labels.csv")

    return annotate_parent, tmp_path / "benchmark_results"


def _run_pipeline(
    annotate_parent: Path,
    benchmark_dir: Path,
    f_measure_top_k: int = 1,
) -> None:
    run_confusion_matrix_step(
        annotate_parent,
        benchmark_dir,
        "truth",
        "cas_cell_type_name_1",
        f_measure_top_k=f_measure_top_k,
    )
    run_aggregate_step(benchmark_dir)
    run_f_measure_step(benchmark_dir)


def test_top_k_flat_f_measure_does_not_change_hierarchical(tmp_path):
    annotate_parent, benchmark_dir = _write_benchmark_fixture(tmp_path)

    _run_pipeline(annotate_parent, benchmark_dir, f_measure_top_k=2)
    run_hierarchical_f_measure_step(benchmark_dir)

    per_sample = pd.read_csv(benchmark_dir / "f_measure_per_sample.csv").iloc[0]
    assert per_sample["tp"] == 3
    assert per_sample["fp"] == 0
    assert per_sample["fn"] == 0
    np.testing.assert_allclose(per_sample["f1_micro"], 1.0)

    per_group = pd.read_csv(benchmark_dir / "f_measure_per_group.csv").iloc[0]
    assert per_group["tp"] == 3
    assert per_group["fp"] == 0
    assert per_group["fn"] == 0
    np.testing.assert_allclose(per_group["f1_micro"], 1.0)

    hierarchical = pd.read_csv(benchmark_dir / "hierarchical_f_measure_per_sample.csv").iloc[0]
    assert hierarchical["h_fp"] > 0
    assert hierarchical["h_fn"] > 0


def test_ranked_inferred_label_columns_requires_top_one_suffix():
    with pytest.raises(ValueError, match="requires --inferred-label to end with '_1'"):
        _ranked_inferred_label_columns("prediction", 2)


def test_top_k_requires_present_ranked_columns(tmp_path):
    annotate_parent, benchmark_dir = _write_benchmark_fixture(tmp_path, include_rank_3=False)

    with pytest.raises(ValueError, match="cas_cell_type_name_3"):
        run_confusion_matrix_step(
            annotate_parent,
            benchmark_dir,
            "truth",
            "cas_cell_type_name_1",
            f_measure_top_k=3,
        )


def test_default_top_one_preserves_existing_flat_f_measure(tmp_path):
    annotate_parent, benchmark_dir = _write_benchmark_fixture(tmp_path)

    _run_pipeline(annotate_parent, benchmark_dir)

    per_sample = pd.read_csv(benchmark_dir / "f_measure_per_sample.csv").iloc[0]
    assert per_sample["tp"] == 1
    assert per_sample["fp"] == 2
    assert per_sample["fn"] == 2
    np.testing.assert_allclose(per_sample["f1_micro"], 1 / 3)

    per_group = pd.read_csv(benchmark_dir / "f_measure_per_group.csv").iloc[0]
    assert per_group["tp"] == 1
    assert per_group["fp"] == 2
    assert per_group["fn"] == 2
    np.testing.assert_allclose(per_group["f1_micro"], 1 / 3)


def test_top_k_option_is_only_on_matrix_build_commands():
    runner = CliRunner()

    confusion_matrix_help = runner.invoke(benchmark_group, ["confusion-matrix", "--help"])
    all_help = runner.invoke(benchmark_group, ["all", "--help"])
    hierarchical_help = runner.invoke(benchmark_group, ["hierarchical", "--help"])

    assert confusion_matrix_help.exit_code == 0
    assert all_help.exit_code == 0
    assert hierarchical_help.exit_code == 0
    assert "--f-measure-top-k" in confusion_matrix_help.output
    assert "--f-measure-top-k" in all_help.output
    assert "--f-measure-top-k" not in hierarchical_help.output
