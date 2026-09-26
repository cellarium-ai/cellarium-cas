"""
Unit tests for azimuth level-spec ordering in map_azimuth_to_cas_labels.
"""

import pandas as pd
import pytest

from cellarium.cas.benchmarking.azimuth.helpers.map_azimuth_to_cas_labels import infer_level_specs


def _azimuth_df(columns: dict) -> pd.DataFrame:
    """Build a 3-cell Azimuth DataFrame from ``{label_col: [values]}`` plus matching score columns."""
    df = pd.DataFrame({col: values for col, values in columns.items()})
    for label_col in columns:
        df[f"{label_col}.score"] = 0.9
    return df


def test_cortex_like_coverage_picks_subclass_first():
    # Column order coarse->fine: class, cluster, subclass, cross_species_cluster.
    # Crosswalk keys are unsuffixed subclass names (Astro, Vip); cluster labels carry
    # suffixes (Astro_1) and cross_species_cluster values are a mix.
    df = _azimuth_df(
        {
            "predicted.class": ["Non-Neuronal", "GABAergic", "GABAergic"],
            "predicted.cluster": ["Astro L1-6 FGFR3 AQP1", "Vip_2", "Sst_3"],
            "predicted.subclass": ["Astro", "Vip", "Sst"],
            "predicted.cross_species_cluster": ["Astro_1", "Vip_2", "Sst_3"],
        }
    )
    crosswalk = {"Astro", "Vip", "Sst", "Astro_1"}
    specs = infer_level_specs(df, crosswalk_labels=crosswalk)
    assert specs[0][0] == "predicted.subclass"
    assert [label for label, _ in specs] == [
        "predicted.subclass",
        "predicted.cross_species_cluster",
        "predicted.cluster",
        "predicted.class",
    ]


def test_bonemarrow_like_coverage_picks_fine_level_first():
    # Column order fine->coarse (l2 then l1); l1 contains 'other' which is not a
    # crosswalk key, while l2 values all map.
    df = _azimuth_df(
        {
            "predicted.celltype.l2": ["Plasma", "Plasma", "CD8 Naive"],
            "predicted.celltype.l1": ["Plasma", "Plasma", "other"],
        }
    )
    crosswalk = {"Plasma", "CD8 Naive"}
    specs = infer_level_specs(df, crosswalk_labels=crosswalk)
    assert [label for label, _ in specs] == ["predicted.celltype.l2", "predicted.celltype.l1"]


def test_coarse_to_fine_tie_breaks_to_finer_level():
    # All levels map fully; coverage ties, so the finer level (more distinct labels)
    # must rank first — same result the old reversal produced for kidney/pbmc/etc.
    df = _azimuth_df(
        {
            "predicted.class": ["Neuron", "Neuron", "Glia"],
            "predicted.subclass": ["Astro", "Vip", "Sst"],
        }
    )
    crosswalk = {"Neuron", "Glia", "Astro", "Vip", "Sst"}
    specs = infer_level_specs(df, crosswalk_labels=crosswalk)
    assert [label for label, _ in specs] == ["predicted.subclass", "predicted.class"]


def test_no_crosswalk_falls_back_to_reversal():
    df = _azimuth_df(
        {
            "predicted.class": ["Neuron", "Neuron", "Glia"],
            "predicted.subclass": ["Astro", "Vip", "Sst"],
        }
    )
    specs = infer_level_specs(df, crosswalk_labels=None)
    assert [label for label, _ in specs] == ["predicted.subclass", "predicted.class"]


def test_zero_coverage_falls_back_to_reversal():
    # No level matches the crosswalk at all; keep current (reversal) behavior.
    df = _azimuth_df(
        {
            "predicted.class": ["Neuron", "Neuron", "Glia"],
            "predicted.subclass": ["Astro", "Vip", "Sst"],
        }
    )
    specs = infer_level_specs(df, crosswalk_labels={"unrelated"})
    assert [label for label, _ in specs] == ["predicted.subclass", "predicted.class"]


def test_raises_without_prediction_pairs():
    df = pd.DataFrame({"foo": [1, 2], "bar": [3, 4]})
    with pytest.raises(ValueError, match="No Azimuth prediction columns"):
        infer_level_specs(df)