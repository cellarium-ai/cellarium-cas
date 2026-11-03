"""
Unit tests for scTab gene-space alignment in the container entrypoint.
"""

import numpy as np
import pytest

from cellarium.cas.benchmarking.comparators.sctab.run import align_to_model_genes, strip_version


def test_strip_version_removes_trailing_digit_suffix():
    assert strip_version("ENSG00000141510.17") == "ENSG00000141510"
    assert strip_version("ENSG1") == "ENSG1"
    assert strip_version("GENE.1.2") == "GENE.1"  # only a trailing .<digits> is stripped


def test_version_strip_matches_and_orders_by_model():
    counts = np.array([[1.0, 2.0], [3.0, 4.0]])
    matrix, matched, expected = align_to_model_genes(
        counts,
        ["ENSG1.3", "ENSG2.7"],
        ["ENSG2", "ENSG1"],
        min_overlap=0.5,
    )
    assert (matched, expected) == (2, 2)
    # Column order follows the model gene list, not the input order.
    np.testing.assert_array_equal(matrix[:, 0], [2.0, 4.0])
    np.testing.assert_array_equal(matrix[:, 1], [1.0, 3.0])


def test_zero_fills_missing_model_genes():
    counts = np.array([[1.0, 2.0], [3.0, 4.0]])
    matrix, matched, expected = align_to_model_genes(
        counts,
        ["A", "C"],
        ["A", "B", "C"],
        min_overlap=0.5,
    )
    assert (matched, expected) == (2, 3)
    np.testing.assert_array_equal(matrix[:, 1], [0.0, 0.0])


def test_duplicate_input_ids_sum_into_one_column_with_single_warning():
    counts = np.array([[1.0, 2.0], [3.0, 4.0]])
    with pytest.warns(UserWarning, match="Duplicate input gene id 'A'"):
        matrix, matched, expected = align_to_model_genes(
            counts,
            ["A.1", "A.2"],
            ["A"],
            min_overlap=0.5,
        )
    assert (matched, expected) == (1, 1)
    np.testing.assert_array_equal(matrix[:, 0], [3.0, 7.0])


def test_below_threshold_raises():
    counts = np.array([[1.0], [3.0]])
    with pytest.raises(ValueError, match="below --min-gene-overlap"):
        align_to_model_genes(
            counts,
            ["ENSG1"],
            ["ENSG1", "ENSG2", "ENSG3", "ENSG4"],
            min_overlap=0.5,
        )
