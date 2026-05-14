import json
import typing as t
import warnings
from pathlib import Path


def map_azimuth_to_cas_labels(
    azimuth_csv_path: str,
    h5ad_path: str,
    output_dir: str,
    crosswalk_csv_path: str,
    crosswalk_azimuth_col: str,
    crosswalk_cl_id_col: str,
    level_specs: t.List[t.Tuple[str, str]],
    azimuth_ref_name: str,
    crosswalk_cl_name_col: t.Optional[str] = None,
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
    :param crosswalk_cl_name_col: Optional column in the crosswalk containing the human-readable
        CL term name (e.g. ``"cl_label"``).  Written to ``cas_cell_type_name_k`` columns.
        If ``None``, the Azimuth label string is used as the cell type name.

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
    obs_names = adata.obs_names.astype(str).tolist()

    azimuth_df = pd.read_csv(azimuth_csv_path, index_col=0)
    azimuth_df.index = azimuth_df.index.astype(str)

    missing_barcodes = set(obs_names) - set(azimuth_df.index)
    if missing_barcodes:
        sample = sorted(missing_barcodes)[:5]
        raise ValueError(
            f"{len(missing_barcodes)} barcodes from h5ad not found in Azimuth CSV " f"(first {len(sample)}): {sample}"
        )
    azimuth_df = azimuth_df.reindex(obs_names)

    missing_level_cols = [
        col for label_col, score_col in level_specs for col in (label_col, score_col) if col not in azimuth_df.columns
    ]
    if missing_level_cols:
        raise ValueError(
            f"The following level columns are not present in the Azimuth CSV: {missing_level_cols}. "
            f"Available columns: {list(azimuth_df.columns)}"
        )

    crosswalk_df = pd.read_csv(crosswalk_csv_path)
    for col in (crosswalk_azimuth_col, crosswalk_cl_id_col):
        if col not in crosswalk_df.columns:
            raise ValueError(
                f"Column '{col}' not found in crosswalk CSV. " f"Available columns: {list(crosswalk_df.columns)}"
            )
    crosswalk_map: t.Dict[str, str] = dict(
        zip(crosswalk_df[crosswalk_azimuth_col].astype(str), crosswalk_df[crosswalk_cl_id_col].astype(str))
    )
    if crosswalk_cl_name_col is not None:
        if crosswalk_cl_name_col not in crosswalk_df.columns:
            raise ValueError(
                f"Column '{crosswalk_cl_name_col}' not found in crosswalk CSV. "
                f"Available columns: {list(crosswalk_df.columns)}"
            )
        crosswalk_name_map: t.Dict[str, str] = dict(
            zip(crosswalk_df[crosswalk_azimuth_col].astype(str), crosswalk_df[crosswalk_cl_name_col].astype(str))
        )
    else:
        crosswalk_name_map = {}

    # --- build label rows ---
    warned_missing: t.Set[str] = set()
    top_k = len(level_specs)
    label_cols = [f"cas_cell_type_label_{k}" for k in range(1, top_k + 1)]
    name_cols = [f"cas_cell_type_name_{k}" for k in range(1, top_k + 1)]
    score_cols = [f"cas_cell_type_score_{k}" for k in range(1, top_k + 1)]
    rows: t.List[t.Dict[str, t.Any]] = []

    for barcode in obs_names:
        row: t.Dict[str, t.Any] = {}
        for rank, (label_col, score_col) in enumerate(level_specs, start=1):
            azimuth_label = azimuth_df.at[barcode, label_col] if label_col in azimuth_df.columns else None
            azimuth_score = azimuth_df.at[barcode, score_col] if score_col in azimuth_df.columns else None

            cl_id = None
            cl_name = None
            if azimuth_label is not None and not pd.isna(azimuth_label):
                label_str = str(azimuth_label)
                cl_id = crosswalk_map.get(label_str)
                if cl_id is None:
                    if label_str not in warned_missing:
                        warnings.warn(
                            f"Azimuth label '{label_str}' (column '{label_col}') not found in crosswalk "
                            f"— filling with 'unknown' for rank {rank}."
                        )
                        warned_missing.add(label_str)
                    cl_id = "unknown"
                    cl_name = "unknown"
                else:
                    cl_name = crosswalk_name_map.get(label_str, label_str)

            score = float(azimuth_score) if azimuth_score is not None and not pd.isna(azimuth_score) else None

            row[f"cas_cell_type_label_{rank}"] = cl_id
            row[f"cas_cell_type_name_{rank}"] = cl_name
            row[f"cas_cell_type_score_{rank}"] = score

        rows.append(row)

    # Column order matches CAS annotate output: label_1..k, name_1..k, score_1..k
    ordered_cols = label_cols + name_cols + score_cols
    labels_df = pd.DataFrame(rows, index=pd.Index(obs_names, name=None))[ordered_cols]
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
