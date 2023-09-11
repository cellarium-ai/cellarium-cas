"""
Test data validation and sanitizing functions.
"""

import typing as t
import unittest

import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp

from cellarium.cas import data_preparation, exceptions

np_random_state = np.random.RandomState(0)


class TestdataPreparation(unittest.TestCase):
    SMALLER_feature_SCHEMA_LENGTH = 12
    INPUT_DATASET_LENGTH = 20

    @staticmethod
    def get_feature_schema() -> t.List:
        return [
            "ENSG00000162639",
            "ENSG00000163568",
            "ENSG00000243701",
            "ENSG00000152580",
            "ENSG00000253811",
            "ENSG00000186439",
            "ENSG00000254756",
            "ENSG00000139211",
            "ENSG00000185615",
            "ENSG00000246379",
            "ENSG00000266469",
            "ENSG00000287710",
            "ENSG00000262006",
        ]

    @staticmethod
    def get_smaller_feature_schema() -> t.List:
        return [
            "ENSG00000162639",
            "ENSG00000163568",
            "ENSG00000243701",
            "ENSG00000253811",
            "ENSG00000186439",
            "ENSG00000254756",
            "ENSG00000139211",
            "ENSG00000185615",
            "ENSG00000246379",
        ]

    @staticmethod
    def get_larger_feature_schema() -> t.List:
        return [
            "ENSG00000162639",
            "ENSG00000163568",
            "ENSG00000243701",
            "ENSG00000152580",
            "ENSG00000253811",
            "ENSG00000186439",
            "ENSG00000254756",
            "ENSG00000139211",
            "ENSG00000185615",
            "ENSG00000263674",
            "ENSG00000246379",
            "ENSG00000266469",
            "ENSG00000287710",
            "ENSG00000185869",
            "ENSG00000262006",
            "ENSG00000167113",
        ]

    @staticmethod
    def get_schema_that_has_missing_and_extra() -> t.List:
        return [
            "ENSG00000163568",
            "ENSG00000243701",
            "ENSG00000152580",
            "ENSG00000253811",
            "ENSG00000186439",
            "ENSG00000254756",
            "ENSG00000139211",
            "ENSG00000263674",
            "ENSG00000246379",
            "ENSG00000266469",
            "ENSG00000287710",
            "ENSG00000167113",
        ]

    @staticmethod
    def get_same_schema_different_order() -> t.List:
        return [
            "ENSG00000162639",
            "ENSG00000243701",
            "ENSG00000163568",
            "ENSG00000152580",
            "ENSG00000186439",
            "ENSG00000253811",
            "ENSG00000254756",
            "ENSG00000139211",
            "ENSG00000246379",
            "ENSG00000185615",
            "ENSG00000266469",
            "ENSG00000262006",
            "ENSG00000287710",
        ]

    def test_data_validation_input_has_missing_features(self):
        cas_feature_schema_list = self.get_feature_schema()
        adata_feature_schema_list = self.get_smaller_feature_schema()
        data_validation_error_was_raised = False
        number_of_missing_features = 0
        number_of_extra_features = 0
        d = np_random_state.randint(0, 500, size=(self.INPUT_DATASET_LENGTH, len(adata_feature_schema_list)))
        adata = anndata.AnnData(
            X=sp.csr_matrix(d),
            obs=pd.DataFrame(index=np.arange(0, self.INPUT_DATASET_LENGTH)),
            var=pd.DataFrame(index=adata_feature_schema_list),
        )
        try:
            data_preparation.validate(
                adata=adata,
                cas_feature_schema_list=cas_feature_schema_list,
                feature_ids_column_name="index",
            )
        except exceptions.DataValidationError as e:
            number_of_missing_features = e.missing_features
            number_of_extra_features = e.extra_features
            data_validation_error_was_raised = True

        self.assertTrue(
            data_validation_error_was_raised,
            msg="`DataValidationError should be raised`",
        )
        self.assertEqual(
            number_of_missing_features,
            len(cas_feature_schema_list) - len(adata_feature_schema_list),
            msg=(
                "Number of missing features should be difference between "
                "`len(cas_feature_schema_list)` and `len(adata_feature_schema_list)`"
            ),
        )
        self.assertEqual(number_of_extra_features, 0, msg="There are no extra features in this test")

    def test_data_validation_input_has_extra_features(self):
        cas_feature_schema_list = self.get_feature_schema()
        adata_feature_schema_list = self.get_larger_feature_schema()
        data_validation_error_was_raised = False
        number_of_missing_features = 0
        number_of_extra_features = 0
        d = np_random_state.randint(0, 500, size=(self.INPUT_DATASET_LENGTH, len(adata_feature_schema_list)))
        adata = anndata.AnnData(
            X=sp.csr_matrix(d),
            obs=pd.DataFrame(index=np.arange(0, self.INPUT_DATASET_LENGTH)),
            var=pd.DataFrame(index=adata_feature_schema_list),
        )
        try:
            data_preparation.validate(
                adata=adata,
                cas_feature_schema_list=cas_feature_schema_list,
                feature_ids_column_name="index",
            )
        except exceptions.DataValidationError as e:
            number_of_missing_features = e.missing_features
            number_of_extra_features = e.extra_features
            data_validation_error_was_raised = True

        self.assertTrue(
            data_validation_error_was_raised,
            msg="`DataValidationError should be raised`",
        )
        self.assertEqual(
            number_of_missing_features,
            0,
            msg="There are no missing values in this test",
        )
        self.assertEqual(
            number_of_extra_features,
            len(adata_feature_schema_list) - len(cas_feature_schema_list),
            msg=(
                "Number of missing features should be difference between "
                "`len(adata_feature_schema_list)` and `len(cas_feature_schema_list)`"
            ),
        )

    def test_data_validation_has_both_missing_and_extra_features(self):
        cas_feature_schema_list = self.get_feature_schema()
        adata_feature_schema_list = self.get_schema_that_has_missing_and_extra()
        data_validation_error_was_raised = False
        number_of_missing_features = 0
        number_of_extra_features = 0
        d = np_random_state.randint(0, 500, size=(self.INPUT_DATASET_LENGTH, len(adata_feature_schema_list)))
        adata = anndata.AnnData(
            X=sp.csr_matrix(d),
            obs=pd.DataFrame(index=np.arange(0, self.INPUT_DATASET_LENGTH)),
            var=pd.DataFrame(index=adata_feature_schema_list),
        )
        try:
            data_preparation.validate(
                adata=adata,
                cas_feature_schema_list=cas_feature_schema_list,
                feature_ids_column_name="index",
            )
        except exceptions.DataValidationError as e:
            number_of_missing_features = e.missing_features
            number_of_extra_features = e.extra_features
            data_validation_error_was_raised = True

        missing_features = len(set(cas_feature_schema_list) - set(adata_feature_schema_list))
        extra_features = len(set(adata_feature_schema_list) - set(cas_feature_schema_list))
        self.assertTrue(
            data_validation_error_was_raised,
            msg="`DataValidationError should be raised`",
        )
        self.assertEqual(number_of_missing_features, missing_features)
        self.assertEqual(number_of_extra_features, extra_features)

    def test_schema_same_but_different_order(self):
        cas_feature_schema_list = self.get_feature_schema()
        adata_feature_schema_list = self.get_same_schema_different_order()
        data_validation_error_was_raised = False
        number_of_missing_features = 0
        number_of_extra_features = 0
        d = np_random_state.randint(0, 500, size=(self.INPUT_DATASET_LENGTH, len(adata_feature_schema_list)))
        adata = anndata.AnnData(
            X=sp.csr_matrix(d),
            obs=pd.DataFrame(index=np.arange(0, self.INPUT_DATASET_LENGTH)),
            var=pd.DataFrame(index=adata_feature_schema_list),
        )
        try:
            data_preparation.validate(
                adata=adata,
                cas_feature_schema_list=cas_feature_schema_list,
                feature_ids_column_name="index",
            )
        except exceptions.DataValidationError as e:
            number_of_missing_features = e.missing_features
            number_of_extra_features = e.extra_features
            data_validation_error_was_raised = True

        self.assertTrue(
            data_validation_error_was_raised,
            msg="`DataValidationError should be raised`",
        )
        self.assertEqual(number_of_missing_features, 0)
        self.assertEqual(number_of_extra_features, 0)

    def test_data_sanitizing(self):
        cas_feature_schema_list = self.get_feature_schema()
        adata_feature_schema_list = self.get_schema_that_has_missing_and_extra()
        d = np_random_state.randint(0, 500, size=(self.INPUT_DATASET_LENGTH, len(adata_feature_schema_list)))
        adata = anndata.AnnData(
            X=sp.csr_matrix(d),
            obs=pd.DataFrame(index=np.arange(0, self.INPUT_DATASET_LENGTH)),
            var=pd.DataFrame(index=adata_feature_schema_list),
        )
        adata_new = data_preparation.sanitize(
            adata=adata,
            cas_feature_schema_list=cas_feature_schema_list,
            count_matrix_name="X",
            feature_ids_column_name="index",
        )
        intersect_features = list(set(cas_feature_schema_list).intersection(set(adata_feature_schema_list)))
        new_matrix_didnt_lose_values = (adata_new[:, intersect_features].X != adata[:, intersect_features].X).sum() == 0
        self.assertTrue(new_matrix_didnt_lose_values)


def main():
    """
    Standard entry point that forwards to unit test main.
    """
    unittest.main()


if __name__ == "__main__":
    main()
