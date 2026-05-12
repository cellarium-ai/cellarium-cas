#!/usr/bin/env python3
"""
Convert Azimuth cell type annotations to the CAS-compatible ``inferred_labels.csv``
and ``metadata.json`` format consumed by ``cellarium-cas benchmark flat``.

Usage
-----
::

    python helpers/map_azimuth_to_cas_labels.py \\
        --azimuth-csv      /data/azimuth_output.csv \\
        --h5ad-path        /data/my_dataset.h5ad \\
        --output-dir       /data/azimuth_annotate_dir \\
        --crosswalk-csv    crosswalk.csv \\
        --crosswalk-azimuth-col  tool_cell_label \\
        --crosswalk-cl-id-col    cl_id \\
        --azimuth-ref-name pbmcref \\
        --level predicted.celltype.l3:predicted.celltype.l3.score \\
        --level predicted.celltype.l2:predicted.celltype.l2.score \\
        --level predicted.celltype.l1:predicted.celltype.l1.score

``--level`` arguments must be ordered **most granular first**.  Each value is
``<azimuth_label_column>:<azimuth_score_column>`` as they appear in the Azimuth
output CSV.

Outputs
-------
- ``<output_dir>/inferred_labels.csv`` — one row per cell, columns
  ``cas_cell_type_label_k`` and ``cas_cell_type_score_k`` for each level k.
- ``<output_dir>/metadata.json`` — provenance info (input path, model name, cell count).
"""

import json
import typing as t
import warnings
from pathlib import Path

import click


def map_azimuth_to_cas_labels(
    azimuth_csv_path: str,
    h5ad_path: str,
    output_dir: str,
    crosswalk_csv_path: str,
    crosswalk_azimuth_col: str,
    crosswalk_cl_id_col: str,
    level_specs: t.List[t.Tuple[str, str]],
    azimuth_ref_name: str,
) -> t.Dict[str, str]:
    """
    Convert Azimuth annotations to CAS-compatible ``inferred_labels.csv`` and ``metadata.json``.

    :param azimuth_csv_path: Path to the Azimuth metadata CSV (row index = barcodes).
    :param h5ad_path: Path to the source ``.h5ad`` file; used for authoritative cell ordering.
    :param output_dir: Directory to write output files into (created if absent).
    :param crosswalk_csv_path: Path to the HRA crosswalk CSV mapping Azimuth labels to CL IDs.
    :param crosswalk_azimuth_col: Column in the crosswalk containing Azimuth cell type labels.
    :param crosswalk_cl_id_col: Column in the crosswalk containing CL ontology term IDs.
    :param level_specs: List of ``(azimuth_label_col, azimuth_score_col)`` tuples ordered
        **most granular first**.  Each tuple names the columns in the Azimuth CSV for one
        annotation level.
    :param azimuth_ref_name: Azimuth reference name (e.g. ``"pbmcref"``), used to build
        the ``model_name`` field as ``azimuth_<ref_name>``.

    :returns: Dict with ``inferred_labels_path`` and ``metadata_path``.
    :raises ValueError: If barcodes in the h5ad are missing from the Azimuth CSV, or if
        required columns are absent from either input CSV.
    """
    import anndata
    import pandas as pd

    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    # --- load inputs ---
    adata = anndata.read_h5ad(h5ad_path)
    obs_names = list(adata.obs_names)

    azimuth_df = pd.read_csv(azimuth_csv_path, index_col=0)
    missing_barcodes = set(obs_names) - set(azimuth_df.index)
    if missing_barcodes:
        sample = sorted(missing_barcodes)[:5]
        raise ValueError(
            f"{len(missing_barcodes)} barcodes from h5ad not found in Azimuth CSV "
            f"(first {len(sample)}): {sample}"
        )
    azimuth_df = azimuth_df.reindex(obs_names)

    crosswalk_df = pd.read_csv(crosswalk_csv_path)
    for col in (crosswalk_azimuth_col, crosswalk_cl_id_col):
        if col not in crosswalk_df.columns:
            raise ValueError(
                f"Column '{col}' not found in crosswalk CSV. "
                f"Available columns: {list(crosswalk_df.columns)}"
            )
    crosswalk_map: t.Dict[str, str] = dict(
        zip(crosswalk_df[crosswalk_azimuth_col].astype(str), crosswalk_df[crosswalk_cl_id_col].astype(str))
    )

    # --- build label rows ---
    warned_missing: t.Set[str] = set()
    rows: t.List[t.Dict[str, t.Any]] = []

    for barcode in obs_names:
        row: t.Dict[str, t.Any] = {}
        for rank, (label_col, score_col) in enumerate(level_specs, start=1):
            azimuth_label = azimuth_df.at[barcode, label_col] if label_col in azimuth_df.columns else None
            azimuth_score = azimuth_df.at[barcode, score_col] if score_col in azimuth_df.columns else None

            cl_id = None
            if azimuth_label is not None and not (
                isinstance(azimuth_label, float) and azimuth_label != azimuth_label  # NaN check
            ):
                label_str = str(azimuth_label)
                cl_id = crosswalk_map.get(label_str)
                if cl_id is None and label_str not in warned_missing:
                    warnings.warn(
                        f"Azimuth label '{label_str}' (column '{label_col}') not found in crosswalk "
                        f"— cell type will be missing for rank {rank}."
                    )
                    warned_missing.add(label_str)

            score = float(azimuth_score) if azimuth_score is not None and azimuth_score == azimuth_score else None

            row[f"cas_cell_type_label_{rank}"] = cl_id
            row[f"cas_cell_type_score_{rank}"] = score

        rows.append(row)

    labels_df = pd.DataFrame(rows, index=pd.Index(obs_names, name="cell_id"))
    labels_path = output_dir_path / "inferred_labels.csv"
    labels_df.to_csv(labels_path)

    # --- metadata ---
    metadata = {
        "input_path": str(Path(h5ad_path).resolve()),
        "model_name": f"azimuth_{azimuth_ref_name}",
        "n_cells": len(obs_names),
        "ontology_resource_name": None,
    }
    metadata_path = output_dir_path / "metadata.json"
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2)

    return {
        "inferred_labels_path": str(labels_path),
        "metadata_path": str(metadata_path),
    }


def _parse_level(value: str) -> t.Tuple[str, str]:
    parts = value.split(":", 1)
    if len(parts) != 2:
        raise click.BadParameter(f"must be in the form LABEL_COL:SCORE_COL, got: '{value}'")
    return parts[0].strip(), parts[1].strip()


@click.command("map-azimuth-to-cas-labels")
@click.option("--azimuth-csv", required=True, type=click.Path(exists=True), help="Path to Azimuth metadata CSV.")
@click.option("--h5ad-path", required=True, type=click.Path(exists=True), help="Path to source .h5ad file.")
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(file_okay=False),
    help="Directory to write output files into.",
)
@click.option("--crosswalk-csv", required=True, type=click.Path(exists=True), help="Path to HRA crosswalk CSV.")
@click.option(
    "--crosswalk-azimuth-col",
    required=True,
    help="Column in crosswalk containing Azimuth cell type labels.",
)
@click.option(
    "--crosswalk-cl-id-col",
    required=True,
    help="Column in crosswalk containing CL ontology term IDs.",
)
@click.option(
    "--azimuth-ref-name",
    required=True,
    help="Azimuth reference name (e.g. pbmcref). Written into model_name as azimuth_<ref_name>.",
)
@click.option(
    "--level",
    "levels",
    multiple=True,
    required=True,
    metavar="LABEL_COL:SCORE_COL",
    help=(
        "Azimuth label and score column pair for one annotation level, as "
        "LABEL_COL:SCORE_COL. Repeat for each level, most granular first."
    ),
)
def main(
    azimuth_csv: str,
    h5ad_path: str,
    output_dir: str,
    crosswalk_csv: str,
    crosswalk_azimuth_col: str,
    crosswalk_cl_id_col: str,
    azimuth_ref_name: str,
    levels: t.Tuple[str, ...],
) -> None:
    """Convert Azimuth annotations to CAS-compatible inferred_labels.csv and metadata.json."""
    level_specs = [_parse_level(lv) for lv in levels]
    try:
        result = map_azimuth_to_cas_labels(
            azimuth_csv_path=azimuth_csv,
            h5ad_path=h5ad_path,
            output_dir=output_dir,
            crosswalk_csv_path=crosswalk_csv,
            crosswalk_azimuth_col=crosswalk_azimuth_col,
            crosswalk_cl_id_col=crosswalk_cl_id_col,
            level_specs=level_specs,
            azimuth_ref_name=azimuth_ref_name,
        )
    except ValueError as e:
        raise click.UsageError(str(e)) from e

    click.echo(f"Saved inferred labels  → {result['inferred_labels_path']}")
    click.echo(f"Saved metadata         → {result['metadata_path']}")


if __name__ == "__main__":
    main()
