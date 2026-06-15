"""
Click command definitions for ``cellarium-cas benchmark``.

Benchmark implementations live in :mod:`cellarium.cas.cli._benchmark_impl`, which imports
optional benchmark dependencies at file level.  This module stays lightweight so docs
tooling can import Click command objects without requiring ``cellarium-cas[benchmark]``.

Pipeline
--------
1. ``cellarium-cas benchmark confusion-matrix``  -- build per-sample CMs  -> cm_raw/
2. ``cellarium-cas benchmark aggregate``          -- aggregate by model    -> cm_aggregate/
3. ``cellarium-cas benchmark f-measure``          -- F-measure CSVs
4. ``cellarium-cas benchmark hierarchical``       -- hierarchical F-measure CSVs
5. ``cellarium-cas benchmark all``                -- convenience: runs all four steps
"""

import click

_annotate_dirs_option = click.option(
    "--annotate-dirs",
    required=True,
    type=click.Path(exists=True),
    help=(
        "Path to a parent directory whose subdirectories are annotate output dirs, "
        "or a .txt file listing one annotate output directory path per line."
    ),
)

_output_dir_option = click.option(
    "--output-dir",
    required=True,
    type=click.Path(file_okay=False),
    help="Benchmarking workspace directory.  All artifacts are written here.",
)

_gt_label_option = click.option(
    "--gt-label",
    required=True,
    help=(
        "Column name in the original .h5ad obs that contains ground-truth cell type labels "
        "(e.g. 'cell_type_ontology_term_id')."
    ),
)

_inferred_label_option = click.option(
    "--inferred-label",
    required=True,
    help=(
        "Column name in inferred_labels.csv that contains the predicted cell type labels "
        "(e.g. 'cas_cell_type_name_1')."
    ),
)

_f_measure_top_k_option = click.option(
    "--f-measure-top-k",
    type=click.IntRange(min=1),
    default=1,
    show_default=True,
    help=(
        "Number of ranked inferred-label columns to consider for flat F-measure. Columns are derived from "
        "--inferred-label by replacing its trailing rank; hierarchical F-measure remains top-1."
    ),
)


@click.group("benchmark")
def benchmark_group() -> None:
    """Compute benchmarking metrics against CAS annotate outputs."""


@benchmark_group.command("confusion-matrix")
@_annotate_dirs_option
@_output_dir_option
@_gt_label_option
@_inferred_label_option
@_f_measure_top_k_option
def confusion_matrix_command(
    annotate_dirs: str,
    output_dir: str,
    gt_label: str,
    inferred_label: str,
    f_measure_top_k: int,
) -> None:
    """Build per-sample confusion matrices and save them to <output-dir>/cm_raw/.

    Reads ground-truth labels from the original .h5ad (via metadata.json input_path)
    and predicted labels from inferred_labels.csv.  The matrix is aligned to the full
    Cell Ontology cl_names universe from ontology_resource.json.
    """
    from ._benchmark_impl import run_confusion_matrix_step

    try:
        result = run_confusion_matrix_step(
            annotate_dirs,
            output_dir,
            gt_label,
            inferred_label,
            f_measure_top_k=f_measure_top_k,
        )
    except (ValueError, FileNotFoundError) as exc:
        raise click.UsageError(str(exc)) from exc

    click.echo(f"Saved {result['n_samples']} confusion matrix(es) -> {result['cm_raw_dir']}")


@benchmark_group.command("aggregate")
@_output_dir_option
def aggregate_command(output_dir: str) -> None:
    """Aggregate per-sample confusion matrices by model name into <output-dir>/cm_aggregate/.

    Groups all matrices in cm_raw/ by their model_name metadata field and sums each group
    into one aggregated sparse confusion matrix.
    """
    from ._benchmark_impl import run_aggregate_step

    try:
        result = run_aggregate_step(output_dir)
    except (ValueError, FileNotFoundError) as exc:
        raise click.UsageError(str(exc)) from exc

    click.echo(f"Aggregated {result['n_groups']} model group(s) -> {result['cm_aggregate_dir']}")


@benchmark_group.command("f-measure")
@_output_dir_option
def f_measure_command(output_dir: str) -> None:
    """Compute standard F-measure metrics from confusion matrices.

    Reads cm_raw/ and cm_aggregate/ and writes:

    \b
      f_measure_per_sample.csv  (columns: model_name, test_sample, tp, fp, fn,
                                           precision_micro, recall_micro, f1_micro, f1_macro)
      f_measure_per_group.csv   (columns: group_name, tp, fp, fn,
                                           precision_micro, recall_micro, f1_micro, f1_macro)
    """
    from ._benchmark_impl import run_f_measure_step

    try:
        result = run_f_measure_step(output_dir)
    except (ValueError, FileNotFoundError) as exc:
        raise click.UsageError(str(exc)) from exc

    click.echo(f"Wrote per-sample F-measure  -> {result['per_sample_path']}")
    click.echo(f"Wrote per-group  F-measure  -> {result['per_group_path']}")


@benchmark_group.command("hierarchical")
@_output_dir_option
def hierarchical_command(output_dir: str) -> None:
    """Compute hierarchical F-measure metrics from confusion matrices.

    Implements the Kiritchenko et al. approach: hierarchical TP/FP/FN are derived
    from the overlap of ontology ancestor sets of the true and predicted labels.
    Reads cm_raw/ and cm_aggregate/ and writes:

    \b
      hierarchical_f_measure_per_sample.csv  (columns: model_name, test_sample,
                                               h_tp, h_fp, h_fn,
                                               h_precision_micro, h_recall_micro,
                                               h_f1_micro, h_f1_macro)
      hierarchical_f_measure_per_group.csv   (columns: group_name, h_tp, h_fp, h_fn,
                                               h_precision_micro, h_recall_micro,
                                               h_f1_micro, h_f1_macro)
    """
    from ._benchmark_impl import run_hierarchical_f_measure_step

    try:
        result = run_hierarchical_f_measure_step(output_dir)
    except (ValueError, FileNotFoundError) as exc:
        raise click.UsageError(str(exc)) from exc

    click.echo(f"Wrote per-sample hierarchical F-measure -> {result['per_sample_path']}")
    click.echo(f"Wrote per-group  hierarchical F-measure -> {result['per_group_path']}")


@benchmark_group.command("all")
@_annotate_dirs_option
@_output_dir_option
@_gt_label_option
@_inferred_label_option
@_f_measure_top_k_option
def all_command(
    annotate_dirs: str,
    output_dir: str,
    gt_label: str,
    inferred_label: str,
    f_measure_top_k: int,
) -> None:
    """Run the full benchmark pipeline in one command.

    Equivalent to running confusion-matrix -> aggregate -> f-measure -> hierarchical
    in sequence.  All artifacts are written to <output-dir>.
    """
    from ._benchmark_impl import run_all_steps

    try:
        result = run_all_steps(
            annotate_dirs,
            output_dir,
            gt_label,
            inferred_label,
            f_measure_top_k=f_measure_top_k,
        )
    except (ValueError, FileNotFoundError) as exc:
        raise click.UsageError(str(exc)) from exc

    click.echo(f"cm_raw/       <- {result['n_samples']} per-sample confusion matrix(es)")
    click.echo(f"cm_aggregate/ <- {result['n_groups']} aggregated model group(s)")
    click.echo(f"F-measure per sample  -> {result['f_measure_per_sample_path']}")
    click.echo(f"F-measure per group   -> {result['f_measure_per_group_path']}")
    click.echo(f"Hierarchical per sample -> {result['h_f_measure_per_sample_path']}")
    click.echo(f"Hierarchical per group  -> {result['h_f_measure_per_group_path']}")
