"""Shared crosswalk-loading helpers for comparator label mapping."""

import typing as t

import pandas as pd


def load_crosswalk_maps(
    crosswalk_csv_path: str,
    crosswalk_label_col: str,
    crosswalk_cl_id_col: str,
    crosswalk_cl_label_col: t.Optional[str],
) -> t.Tuple[t.Dict[str, str], t.Dict[str, str]]:
    """Load crosswalk CSV and return ``(crosswalk_map, crosswalk_label_map)``."""
    crosswalk_df = pd.read_csv(crosswalk_csv_path)
    required_cols = [crosswalk_label_col, crosswalk_cl_id_col]
    if crosswalk_cl_label_col is not None:
        required_cols.append(crosswalk_cl_label_col)
    for col in required_cols:
        if col not in crosswalk_df.columns:
            raise ValueError(
                f"Column '{col}' not found in crosswalk CSV. Available columns: {list(crosswalk_df.columns)}"
            )
    crosswalk_map: t.Dict[str, str] = dict(
        zip(crosswalk_df[crosswalk_label_col].astype(str), crosswalk_df[crosswalk_cl_id_col].astype(str))
    )
    if crosswalk_cl_label_col is not None:
        crosswalk_label_map: t.Dict[str, str] = dict(
            zip(crosswalk_df[crosswalk_label_col].astype(str), crosswalk_df[crosswalk_cl_label_col].astype(str))
        )
    else:
        crosswalk_label_map = {}
    return crosswalk_map, crosswalk_label_map
