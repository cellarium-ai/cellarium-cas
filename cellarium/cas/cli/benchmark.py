"""
Click command definitions for ``cellarium-cas benchmark``.

Benchmark implementations live in :mod:`cellarium.cas.cli._benchmark_impl`, which imports
optional benchmark dependencies at file level. This module stays lightweight so docs tooling
can import Click command objects without requiring ``cellarium-cas[benchmark]``.
"""

import typing as t
from pathlib import Path

import click


def run_flat_benchmark(
    annotate_dirs: t.Union[str, Path],
    gt_column_name: str,
    output_dir: t.Union[str, Path],
    save_class_level: bool = False,
) -> t.Dict[str, t.Any]:
    """Compute flat classification metrics across one or more annotate output directories."""
    from ._benchmark_impl import run_flat_benchmark as _run_flat_benchmark

    return _run_flat_benchmark(
        annotate_dirs=annotate_dirs,
        gt_column_name=gt_column_name,
        output_dir=output_dir,
        save_class_level=save_class_level,
    )


def run_ontology_aware_benchmark(
    annotate_dirs: t.Union[str, Path],
    gt_cl_column_name: str,
    output_dir: t.Union[str, Path],
    num_hops: int = 4,
    save_cell_level: bool = False,
) -> t.Dict[str, t.Any]:
    """Compute ontology-aware metrics across one or more annotate output directories."""
    from ._benchmark_impl import run_ontology_aware_benchmark as _run_ontology_aware_benchmark

    return _run_ontology_aware_benchmark(
        annotate_dirs=annotate_dirs,
        gt_cl_column_name=gt_cl_column_name,
        output_dir=output_dir,
        num_hops=num_hops,
        save_cell_level=save_cell_level,
    )


def run_hierarchical_f_measure_benchmark(
    annotate_dirs: t.Union[str, Path],
    gt_cl_column_name: str,
    output_dir: t.Union[str, Path],
    save_class_level: bool = False,
) -> t.Dict[str, t.Any]:
    """Compute hierarchical F-measure metrics across one or more annotate output directories."""
    from ._benchmark_impl import run_hierarchical_f_measure_benchmark as _run_hierarchical_f_measure_benchmark

    return _run_hierarchical_f_measure_benchmark(
        annotate_dirs=annotate_dirs,
        gt_cl_column_name=gt_cl_column_name,
        output_dir=output_dir,
        save_class_level=save_class_level,
    )


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
    "--save-class-level/--no-save-class-level",
    default=False,
    show_default=True,
    help="Also save one class-level flat metric CSV with sample and total rows.",
)
def flat_command(
    annotate_dirs: str,
    gt_column_name: str,
    output_dir: str,
    save_class_level: bool,
) -> None:
    """Compute flat classification metrics across one or more annotate output directories."""
    try:
        result = run_flat_benchmark(
            annotate_dirs=annotate_dirs,
            gt_column_name=gt_column_name,
            output_dir=output_dir,
            save_class_level=save_class_level,
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
