"""
CLI subcommand group: cellarium-cas benchmark

Subcommands:
  flat             -- Compute flat (taxonomy-agnostic) classification metrics.
  ontology-aware   -- Compute ontology-aware metrics using hop-distance scoring.

The core logic lives in :func:`run_flat_benchmark` and :func:`run_ontology_aware_benchmark`,
which can be imported and called directly from Python (e.g. KFP components).
The click commands are thin wrappers around them.
"""

import re
import typing as t
from pathlib import Path

import click
import pandas as pd

from cellarium.cas.benchmarking.flat import compute_flat_metrics
from cellarium.cas.benchmarking.hierarchical_f_measure import compute_hierarchical_f_measure_metrics
from cellarium.cas.benchmarking.ontology_aware import compute_ontology_aware_metrics
from cellarium.cas.models import CellTypeOntologyAwareResults

from ._io import (
    METADATA_FILENAME,
    ONTOLOGY_RESOURCE_FILENAME,
    ONTOLOGY_RESPONSE_FILENAME,
    collect_annotate_output_dirs,
    load_inferred_labels,
    load_metadata,
    load_ontology_resource,
    load_ontology_response,
)


def run_flat_benchmark(
    annotate_dirs: t.Union[str, Path],
    gt_column_name: str,
    output_dir: t.Union[str, Path],
    save_cell_level: bool = False,
) -> t.Dict[str, t.Any]:
    """
    Compute flat classification metrics across one or more annotate output directories.

    This is the plain-Python entry point. The CLI command :func:`flat_command` is a thin
    wrapper around this function.

    :param annotate_dirs: Path to a parent directory whose subdirectories are annotate output
        directories, or a ``.txt`` file listing one annotate output directory path per line.
    :param gt_column_name: ``obs`` column name containing flat ground truth labels in the
        original ``.h5ad`` file (path read from each directory's ``metadata.json``).
    :param output_dir: Directory to write ``flat_benchmark_summary.csv`` into.
    :param save_cell_level: If ``True``, also save per-directory cell-level CSVs.

    :returns: Dict with ``summary_path``, ``n_combinations``, ``n_summary_rows``.
    :raises ValueError: If a required obs column is missing or no valid directories are found.
    """
    import anndata

    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    dirs = collect_annotate_output_dirs(annotate_dirs)
    summary_rows: t.List[pd.DataFrame] = []

    for d in dirs:
        metadata = load_metadata(d)
        input_path = metadata["input_path"]
        model_name = metadata.get("model_name", "unknown")

        adata = anndata.read_h5ad(input_path)
        if gt_column_name not in adata.obs.columns:
            raise ValueError(
                f"Ground truth column '{gt_column_name}' not found in obs of {input_path}. "
                f"Available columns: {list(adata.obs.columns)}"
            )
        ground_truths = list(adata.obs[gt_column_name].astype(str))

        labels_df = load_inferred_labels(d)
        top_k = _infer_top_k(labels_df)
        predictions = _extract_label_predictions(labels_df, top_k)

        summary = compute_flat_metrics(ground_truths, predictions, top_k=top_k, cell_level=False)
        summary.insert(0, "model_name", model_name)
        summary.insert(0, "input_path", input_path)
        summary.insert(0, "annotate_dir", d.name)
        summary_rows.append(summary)

        if save_cell_level:
            cell_df = compute_flat_metrics(ground_truths, predictions, top_k=top_k, cell_level=True)
            cell_path = output_dir_path / f"{d.name}_cell_level_flat.csv"
            cell_df.to_csv(cell_path, index=False)

    if not summary_rows:
        raise ValueError("No results produced — check that annotate_dirs contains valid directories.")

    combined = pd.concat(summary_rows, ignore_index=True)
    summary_path = output_dir_path / "flat_benchmark_summary.csv"
    combined.to_csv(summary_path, index=False)

    return {
        "summary_path": str(summary_path),
        "n_combinations": len(dirs),
        "n_summary_rows": len(combined),
    }


def run_ontology_aware_benchmark(
    annotate_dirs: t.Union[str, Path],
    gt_cl_column_name: str,
    output_dir: t.Union[str, Path],
    num_hops: int = 4,
    save_cell_level: bool = False,
) -> t.Dict[str, t.Any]:
    """
    Compute ontology-aware metrics across one or more annotate output directories.

    This is the plain-Python entry point. The CLI command :func:`ontology_aware_command` is a
    thin wrapper around this function.

    :param annotate_dirs: Path to a parent directory whose subdirectories are annotate output
        directories, or a ``.txt`` file listing one annotate output directory path per line.
        Each directory must contain ``ontology_resource.json`` (saved by ``annotate`` with
        ``save_ontology_resource=True``).
    :param gt_cl_column_name: ``obs`` column name containing Cell Ontology term IDs
        (e.g. ``"CL:0000121"``) in the original ``.h5ad`` file.
    :param output_dir: Directory to write ``ontology_aware_benchmark_summary.csv`` into.
    :param num_hops: Maximum hop distance for ontology-aware scoring.
    :param save_cell_level: If ``True``, also save per-directory cell-level CSVs.

    :returns: Dict with ``summary_path`` and ``n_combinations``.
    :raises ValueError: If a required obs column is missing or ``ontology_resource.json`` is absent.
    """
    import anndata

    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    dirs = collect_annotate_output_dirs(annotate_dirs)
    summary_rows: t.List[pd.DataFrame] = []

    for d in dirs:
        metadata = load_metadata(d)
        input_path = metadata["input_path"]
        model_name = metadata.get("model_name", "unknown")

        adata = anndata.read_h5ad(input_path)
        if gt_cl_column_name not in adata.obs.columns:
            raise ValueError(
                f"Ground truth column '{gt_cl_column_name}' not found in obs of {input_path}. "
                f"Available columns: {list(adata.obs.columns)}"
            )
        ground_truths = list(adata.obs[gt_cl_column_name].astype(str))

        response = load_ontology_response(d)
        ontology_resource = load_ontology_resource(d)

        summary = compute_ontology_aware_metrics(
            response=response,
            ground_truths=ground_truths,
            ontology_resource=ontology_resource,
            num_hops=num_hops,
            cell_level=False,
        )
        summary.insert(0, "model_name", model_name)
        summary.insert(0, "input_path", input_path)
        summary.insert(0, "annotate_dir", d.name)
        summary_rows.append(summary)

        if save_cell_level:
            cell_df = compute_ontology_aware_metrics(
                response=response,
                ground_truths=ground_truths,
                ontology_resource=ontology_resource,
                num_hops=num_hops,
                cell_level=True,
            )
            cell_path = output_dir_path / f"{d.name}_cell_level_ontology_aware.csv"
            cell_df.to_csv(cell_path, index=False)

    if not summary_rows:
        raise ValueError("No results produced — check that annotate_dirs contains valid directories.")

    combined = pd.concat(summary_rows, ignore_index=True)
    summary_path = output_dir_path / "ontology_aware_benchmark_summary.csv"
    combined.to_csv(summary_path, index=False)

    return {
        "summary_path": str(summary_path),
        "n_combinations": len(dirs),
    }


def run_hierarchical_f_measure_benchmark(
    annotate_dirs: t.Union[str, Path],
    gt_cl_column_name: str,
    output_dir: t.Union[str, Path],
    save_class_level: bool = False,
) -> t.Dict[str, t.Any]:
    """
    Compute hierarchical F-measure metrics across one or more annotate output directories.

    This is the plain-Python entry point. The CLI command :func:`hierarchical_f_measure_command`
    is a thin wrapper around this function.
    """
    import anndata

    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    dirs = collect_annotate_output_dirs(
        annotate_dirs,
        required_files=(ONTOLOGY_RESPONSE_FILENAME, ONTOLOGY_RESOURCE_FILENAME, METADATA_FILENAME),
        missing_hint="Ensure 'cellarium-cas annotate' was run with --save-metadata and --save-ontology-resource.",
    )
    summary_rows: t.List[pd.DataFrame] = []
    class_rows: t.List[pd.DataFrame] = []
    totals_by_model: t.Dict[str, t.Dict[str, t.Any]] = {}

    def _add_hierarchical_context_columns(
        df: pd.DataFrame, annotate_dir: str, row_type: str, input_path: str, model_name: str
    ) -> pd.DataFrame:
        df = df.copy()
        df.insert(0, "model_name", model_name)
        df.insert(0, "input_path", input_path)
        df.insert(0, "row_type", row_type)
        df.insert(0, "annotate_dir", annotate_dir)
        return df

    for d in dirs:
        metadata = load_metadata(d)
        input_path = metadata["input_path"]
        model_name = metadata.get("model_name", "unknown")

        adata = anndata.read_h5ad(input_path)
        if gt_cl_column_name not in adata.obs.columns:
            raise ValueError(
                f"Ground truth column '{gt_cl_column_name}' not found in obs of {input_path}. "
                f"Available columns: {list(adata.obs.columns)}"
            )
        ground_truths = list(adata.obs[gt_cl_column_name].astype(str))

        response = load_ontology_response(d)
        ontology_resource = load_ontology_resource(d)

        total = totals_by_model.setdefault(
            model_name,
            {
                "annotations": [],
                "ground_truths": [],
                "ontology_resource": ontology_resource,
            },
        )
        if total["ontology_resource"] != ontology_resource:
            raise ValueError(
                f"All annotate output directories for model '{model_name}' must use the same ontology resource "
                f"to compute global hierarchical F-measure totals."
            )
        total["annotations"].extend(response.data)
        total["ground_truths"].extend(ground_truths)

        summary = compute_hierarchical_f_measure_metrics(
            response=response,
            ground_truths=ground_truths,
            ontology_resource=ontology_resource,
            class_level=False,
        )
        summary_rows.append(_add_hierarchical_context_columns(summary, d.name, "sample", input_path, model_name))

        if save_class_level:
            class_df = compute_hierarchical_f_measure_metrics(
                response=response,
                ground_truths=ground_truths,
                ontology_resource=ontology_resource,
                class_level=True,
            )
            class_rows.append(_add_hierarchical_context_columns(class_df, d.name, "sample", input_path, model_name))

    for model_name in sorted(totals_by_model):
        total = totals_by_model[model_name]
        response = CellTypeOntologyAwareResults(data=total["annotations"], model_name=model_name)
        summary = compute_hierarchical_f_measure_metrics(
            response=response,
            ground_truths=total["ground_truths"],
            ontology_resource=total["ontology_resource"],
            class_level=False,
        )
        summary_rows.append(_add_hierarchical_context_columns(summary, "__total__", "total", "__all__", model_name))

        if save_class_level:
            class_df = compute_hierarchical_f_measure_metrics(
                response=response,
                ground_truths=total["ground_truths"],
                ontology_resource=total["ontology_resource"],
                class_level=True,
            )
            class_rows.append(_add_hierarchical_context_columns(class_df, "__total__", "total", "__all__", model_name))

    if not summary_rows:
        raise ValueError("No results produced — check that annotate_dirs contains valid directories.")

    combined = pd.concat(summary_rows, ignore_index=True)
    summary_path = output_dir_path / "hierarchical_f_measure_summary.csv"
    combined.to_csv(summary_path, index=False)

    class_level_path = None
    n_class_level_rows = 0
    if save_class_level:
        class_level = pd.concat(class_rows, ignore_index=True) if class_rows else pd.DataFrame()
        class_level_path = output_dir_path / "hierarchical_f_measure_class_level.csv"
        class_level.to_csv(class_level_path, index=False)
        n_class_level_rows = len(class_level)

    result = {
        "summary_path": str(summary_path),
        "n_combinations": len(dirs),
        "n_summary_rows": len(combined),
    }
    if class_level_path is not None:
        result["class_level_path"] = str(class_level_path)
        result["n_class_level_rows"] = n_class_level_rows
    return result


def _infer_top_k(labels_df: pd.DataFrame) -> int:
    """Infer top_k by counting ``cas_cell_type_label_k`` columns in *labels_df*."""
    k = 0
    for col in labels_df.columns:
        m = re.match(r"^cas_cell_type_label_(\d+)$", col)
        if m:
            k = max(k, int(m.group(1)))
    if k == 0:
        raise ValueError(
            "Could not find any 'cas_cell_type_label_k' columns in inferred_labels.csv. "
            "Ensure 'cellarium-cas annotate' was run with --infer-labels."
        )
    return k


def _extract_label_predictions(labels_df: pd.DataFrame, top_k: int) -> t.List[t.List[str]]:
    """Return ranked label predictions as a list-of-lists (one per cell) from *labels_df*."""
    predictions: t.List[t.List[str]] = []
    for _, row in labels_df.iterrows():
        cell_preds = []
        for k in range(1, top_k + 1):
            col = f"cas_cell_type_label_{k}"
            if col in labels_df.columns:
                val = row[col]
                if pd.notna(val):
                    cell_preds.append(str(val))
        predictions.append(cell_preds)
    return predictions


# CLI commands (thin wrappers over the plain functions above)


@click.group("benchmark")
def benchmark_group() -> None:
    """Compute benchmarking metrics against CAS annotate outputs."""


@benchmark_group.command("flat")
@click.option(
    "--annotate-dirs",
    required=True,
    type=click.Path(exists=True),
    help=(
        "Path to a parent directory whose subdirectories are annotate output dirs, "
        "or a .txt file listing one annotate output directory path per line."
    ),
)
@click.option(
    "--gt-column-name", required=True, help="obs column name containing flat ground truth labels in the original .h5ad."
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(file_okay=False),
    help="Directory to write benchmark summary CSV into.",
)
@click.option(
    "--save-cell-level/--no-save-cell-level",
    default=False,
    show_default=True,
    help="Also save per-directory cell-level metric CSVs.",
)
def flat_command(
    annotate_dirs: str,
    gt_column_name: str,
    output_dir: str,
    save_cell_level: bool,
) -> None:
    """Compute flat classification metrics across one or more annotate output directories."""
    try:
        result = run_flat_benchmark(
            annotate_dirs=annotate_dirs,
            gt_column_name=gt_column_name,
            output_dir=output_dir,
            save_cell_level=save_cell_level,
        )
    except ValueError as e:
        raise click.UsageError(str(e)) from e

    click.echo(f"Saved flat benchmark summary ({result['n_summary_rows']} rows) → {result['summary_path']}")


@benchmark_group.command("ontology-aware")
@click.option(
    "--annotate-dirs",
    required=True,
    type=click.Path(exists=True),
    help=(
        "Path to a parent directory whose subdirectories are annotate output dirs, "
        "or a .txt file listing one annotate output directory path per line."
    ),
)
@click.option(
    "--gt-cl-column-name",
    required=True,
    help="obs column name containing Cell Ontology term IDs (e.g. CL:0000121) in the original .h5ad.",
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(file_okay=False),
    help="Directory to write benchmark summary CSV into.",
)
@click.option("--num-hops", default=4, show_default=True, help="Maximum hop distance for ontology-aware scoring.")
@click.option(
    "--save-cell-level/--no-save-cell-level",
    default=False,
    show_default=True,
    help="Also save per-directory cell-level metric CSVs.",
)
def ontology_aware_command(
    annotate_dirs: str,
    gt_cl_column_name: str,
    output_dir: str,
    num_hops: int,
    save_cell_level: bool,
) -> None:
    """Compute ontology-aware metrics across one or more annotate output directories."""
    try:
        result = run_ontology_aware_benchmark(
            annotate_dirs=annotate_dirs,
            gt_cl_column_name=gt_cl_column_name,
            output_dir=output_dir,
            num_hops=num_hops,
            save_cell_level=save_cell_level,
        )
    except (ValueError, FileNotFoundError) as e:
        raise click.UsageError(str(e)) from e

    click.echo(f"Saved ontology-aware benchmark summary → {result['summary_path']}")


@benchmark_group.command("hierarchical-f-measure")
@click.option(
    "--annotate-dirs",
    required=True,
    type=click.Path(exists=True),
    help=(
        "Path to a parent directory whose subdirectories are annotate output dirs, "
        "or a .txt file listing one annotate output directory path per line."
    ),
)
@click.option(
    "--gt-cl-column-name",
    required=True,
    help="obs column name containing Cell Ontology term IDs (e.g. CL:0000121) in the original .h5ad.",
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(file_okay=False),
    help="Directory to write benchmark summary CSV into.",
)
@click.option(
    "--save-class-level/--no-save-class-level",
    default=False,
    show_default=True,
    help="Also save one class-level hierarchical F-measure CSV with sample and total rows.",
)
def hierarchical_f_measure_command(
    annotate_dirs: str,
    gt_cl_column_name: str,
    output_dir: str,
    save_class_level: bool,
) -> None:
    """Compute hierarchical F-measure metrics across one or more annotate output directories."""
    try:
        result = run_hierarchical_f_measure_benchmark(
            annotate_dirs=annotate_dirs,
            gt_cl_column_name=gt_cl_column_name,
            output_dir=output_dir,
            save_class_level=save_class_level,
        )
    except (ValueError, FileNotFoundError) as e:
        raise click.UsageError(str(e)) from e

    click.echo(f"Saved hierarchical F-measure summary ({result['n_summary_rows']} rows) → {result['summary_path']}")
