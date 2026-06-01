import json
import typing as t
import warnings
from pathlib import Path


def _load_and_align_azimuth_df(azimuth_csv_path: str, obs_names: t.List[str]) -> t.Any:
    """Load the Azimuth CSV, validate barcode coverage, and reindex to *obs_names* order."""
    import pandas as pd

    azimuth_df = pd.read_csv(azimuth_csv_path, index_col=0)
    azimuth_df.index = azimuth_df.index.astype(str)
    missing_barcodes = set(obs_names) - set(azimuth_df.index)
    if missing_barcodes:
        sample = sorted(missing_barcodes)[:5]
        raise ValueError(
            f"{len(missing_barcodes)} barcodes from h5ad not found in Azimuth CSV (first {len(sample)}): {sample}"
        )
    return azimuth_df.reindex(obs_names)


def _validate_level_cols(level_specs: t.List[t.Tuple[str, str]], azimuth_df: t.Any) -> None:
    """Raise ``ValueError`` if any label/score column from *level_specs* is absent in *azimuth_df*."""
    missing = [
        col for label_col, score_col in level_specs for col in (label_col, score_col) if col not in azimuth_df.columns
    ]
    if missing:
        raise ValueError(
            f"The following level columns are not present in the Azimuth CSV: {missing}. "
            f"Available columns: {list(azimuth_df.columns)}"
        )


def _load_crosswalk_maps(
    crosswalk_csv_path: str,
    crosswalk_azimuth_col: str,
    crosswalk_cl_id_col: str,
    crosswalk_cl_label_col: t.Optional[str],
) -> t.Tuple[t.Dict[str, str], t.Dict[str, str]]:
    """Load crosswalk CSV and return ``(crosswalk_map, crosswalk_label_map)``."""
    import pandas as pd

    crosswalk_df = pd.read_csv(crosswalk_csv_path)
    required_cols = [crosswalk_azimuth_col, crosswalk_cl_id_col]
    if crosswalk_cl_label_col is not None:
        required_cols.append(crosswalk_cl_label_col)
    for col in required_cols:
        if col not in crosswalk_df.columns:
            raise ValueError(
                f"Column '{col}' not found in crosswalk CSV. Available columns: {list(crosswalk_df.columns)}"
            )
    crosswalk_map: t.Dict[str, str] = dict(
        zip(crosswalk_df[crosswalk_azimuth_col].astype(str), crosswalk_df[crosswalk_cl_id_col].astype(str))
    )
    if crosswalk_cl_label_col is not None:
        crosswalk_label_map: t.Dict[str, str] = dict(
            zip(crosswalk_df[crosswalk_azimuth_col].astype(str), crosswalk_df[crosswalk_cl_label_col].astype(str))
        )
    else:
        crosswalk_label_map = {}
    return crosswalk_map, crosswalk_label_map


def _map_azimuth_label(
    azimuth_label: t.Any,
    label_col: str,
    rank: int,
    crosswalk_map: t.Dict[str, str],
    crosswalk_label_map: t.Dict[str, str],
    warned_missing: t.Set[str],
) -> t.Tuple[t.Optional[str], t.Optional[str]]:
    """Map a single Azimuth label to ``(cl_id, cl_label)``, warning once per unknown label."""
    import pandas as pd

    if azimuth_label is None or pd.isna(azimuth_label):
        return None, None
    label_str = str(azimuth_label)
    cl_id = crosswalk_map.get(label_str)
    if cl_id is not None:
        return cl_id, crosswalk_label_map.get(label_str, label_str)
    if label_str not in warned_missing:
        warnings.warn(
            f"Azimuth label '{label_str}' (column '{label_col}') not found in crosswalk "
            f"— filling with 'unknown' for rank {rank}."
        )
        warned_missing.add(label_str)
    return "unknown", "unknown"


def _build_cell_row(
    barcode: str,
    level_specs: t.List[t.Tuple[str, str]],
    azimuth_df: t.Any,
    crosswalk_map: t.Dict[str, str],
    crosswalk_label_map: t.Dict[str, str],
    warned_missing: t.Set[str],
) -> t.Dict[str, t.Any]:
    """Build one output row dict for *barcode* across all annotation levels."""
    import pandas as pd

    row: t.Dict[str, t.Any] = {}
    for rank, (label_col, score_col) in enumerate(level_specs, start=1):
        azimuth_label = azimuth_df.at[barcode, label_col]
        azimuth_score = azimuth_df.at[barcode, score_col]
        cl_id, cl_label = _map_azimuth_label(
            azimuth_label, label_col, rank, crosswalk_map, crosswalk_label_map, warned_missing
        )
        score = float(azimuth_score) if azimuth_score is not None and not pd.isna(azimuth_score) else None
        row[f"cas_cell_type_label_{rank}"] = cl_label
        row[f"cas_cell_type_name_{rank}"] = cl_id
        row[f"cas_cell_type_score_{rank}"] = score
    return row


def infer_level_specs(azimuth_df: t.Any) -> t.List[t.Tuple[str, str]]:
    """
    Auto-detect ``(label_col, score_col)`` pairs from an Azimuth output DataFrame.

    Azimuth writes predicted labels as ``predicted.<level>`` and per-cell confidence
    scores as ``predicted.<level>.score``.  Pairs are returned **coarsest first**
    (column-position order, since Azimuth adds coarser levels first), so rank 1 is always populated.

    :param azimuth_df: DataFrame loaded from the Azimuth metadata CSV.
    :returns: List of ``(label_col, score_col)`` tuples, coarsest first.
    :raises ValueError: If no ``predicted.*`` / ``predicted.*.score`` pairs are found.
    """
    pairs = [
        (col, f"{col}.score")
        for col in azimuth_df.columns
        if col.startswith("predicted.") and not col.endswith(".score") and f"{col}.score" in azimuth_df.columns
    ]
    if not pairs:
        raise ValueError(
            "No Azimuth prediction columns found.  Expected columns matching "
            "'predicted.<level>' with a corresponding 'predicted.<level>.score' column.  "
            f"Available columns: {list(azimuth_df.columns)}"
        )
    # Azimuth adds annotation levels coarse-to-fine; preserve that order so rank 1 is always populated.
    return pairs


def map_azimuth_to_cas_labels(
    azimuth_csv_path: str,
    h5ad_path: str,
    output_dir: str,
    crosswalk_csv_path: str,
    crosswalk_azimuth_col: str,
    crosswalk_cl_id_col: str,
    azimuth_ref_name: str,
    level_specs: t.Optional[t.List[t.Tuple[str, str]]] = None,
    crosswalk_cl_label_col: t.Optional[str] = None,
) -> t.Dict[str, str]:
    """
    Convert Azimuth annotations to CAS-compatible ``inferred_labels.csv`` and ``metadata.json``.

    :param azimuth_csv_path: Path to the Azimuth metadata CSV (row index = barcodes).
    :param h5ad_path: Path to the source ``.h5ad`` file; used for authoritative cell ordering.
    :param output_dir: Directory to write output files into (created if absent).
    :param crosswalk_csv_path: Path to the HRA crosswalk CSV mapping Azimuth labels to CL IDs.
    :param crosswalk_azimuth_col: Column in the crosswalk containing Azimuth cell type labels.
    :param crosswalk_cl_id_col: Column in the crosswalk containing CL ontology term IDs
        (e.g. ``"CL:0000084"``).  Written to ``cas_cell_type_name_k`` columns.
    :param azimuth_ref_name: Azimuth reference name (e.g. ``"pbmcref"``), used to build
        the ``model_name`` field as ``azimuth_<ref_name>``.
    :param level_specs: List of ``(azimuth_label_col, azimuth_score_col)`` tuples ordered
        **coarsest first**.  Each tuple names the columns in the Azimuth CSV for one
        annotation level.  If ``None`` (default), pairs are auto-detected from columns
        matching ``predicted.<level>`` / ``predicted.<level>.score``.
    :param crosswalk_cl_label_col: Optional column in the crosswalk containing the human-readable
        CL term label (e.g. ``"cl_label"``).  Written to ``cas_cell_type_label_k`` columns.
        If ``None``, the Azimuth label string is used.

    :returns: Dict with ``inferred_labels_path`` and ``metadata_path``.
    :raises ValueError: If barcodes in the h5ad are missing from the Azimuth CSV, or if
        required columns are absent from either input CSV.
    """
    import anndata
    import pandas as pd

    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    adata = anndata.read_h5ad(h5ad_path)
    obs_names = adata.obs_names.astype(str).tolist()

    azimuth_df = _load_and_align_azimuth_df(azimuth_csv_path, obs_names)
    if level_specs is None:
        level_specs = infer_level_specs(azimuth_df)
    _validate_level_cols(level_specs, azimuth_df)

    crosswalk_map, crosswalk_label_map = _load_crosswalk_maps(
        crosswalk_csv_path, crosswalk_azimuth_col, crosswalk_cl_id_col, crosswalk_cl_label_col
    )

    # --- build label rows ---
    warned_missing: t.Set[str] = set()
    top_k = len(level_specs)
    label_cols = [f"cas_cell_type_label_{k}" for k in range(1, top_k + 1)]
    name_cols = [f"cas_cell_type_name_{k}" for k in range(1, top_k + 1)]
    score_cols = [f"cas_cell_type_score_{k}" for k in range(1, top_k + 1)]

    rows = [
        _build_cell_row(barcode, level_specs, azimuth_df, crosswalk_map, crosswalk_label_map, warned_missing)
        for barcode in obs_names
    ]

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
