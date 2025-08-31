import os
import unittest
from unittest.mock import patch
from PIL import Image, ImageChops
import numpy as np
import pandas as pd
from module.core.FileSystem import FileSystem
from module.core.plotters import (
    histogram,
    summary_histogram,
    correlogram,
    network,
    network_summary,
    network_degrees,
    correlation,
    statistics_table,
)


class TestPlotters(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.project_name = "tcb2"
        cls.dataset_name = "hplc"

        # Patch where the names are actually used
        cls.p_yes_meta = patch("module.core.Metadata.yes_or_no", return_value=True)
        cls.p_yes_ds = patch("module.core.ProjectDataset.yes_or_no", return_value=True)
        cls.p_input = patch(
            "module.core.ProjectDataset.input_escape",
            side_effect=[
                "module/example_project/tcb2_hplc_data.csv",
                "module/example_project/tcb2_fake_behavior_data.csv",
            ],
        )
        cls.p_edit = patch(
            "module.core.Metadata.ProjectMetadata.make_user_edit_excel",
            return_value=None,
        )

        # Start patches
        cls.m_yes_meta = cls.p_yes_meta.start()
        cls.m_yes_ds = cls.p_yes_ds.start()
        cls.m_input = cls.p_input.start()
        cls.m_edit = cls.p_edit.start()

    @classmethod
    def tearDownClass(cls):
        cls.p_yes_meta.stop()
        cls.p_yes_ds.stop()
        cls.p_input.stop()
        cls.p_edit.stop()
        FileSystem.delete_project(cls.project_name)

    def assert_image_similar(self, expected_path, actual_path, threshold=0):
        """
        Assert that two images are similar within the given threshold.

        Args:
            expected_path: Path to expected image
            actual_path: Path to actual image
            threshold: Percentage difference threshold (0.01 = 1%)
        """
        self.assertTrue(
            os.path.exists(expected_path), f"Expected image not found: {expected_path}"
        )
        self.assertTrue(
            os.path.exists(actual_path), f"Actual image not found: {actual_path}"
        )

        try:
            with (
                Image.open(expected_path) as expected,
                Image.open(actual_path) as actual,
            ):
                # Ensure images are same size
                if expected.size != actual.size:
                    self.fail(
                        f"Image sizes differ: expected {expected.size}, got {actual.size}"
                    )

                # Convert to same mode if needed
                if expected.mode != actual.mode:
                    actual = actual.convert(expected.mode)

                # Calculate difference
                diff = ImageChops.difference(expected, actual)
                diff_array = np.array(diff)

                # Calculate percentage of different pixels
                if len(diff_array.shape) == 3:  # RGB
                    different_pixels = np.any(diff_array > 0, axis=2)
                else:  # Grayscale
                    different_pixels = diff_array > 0

                diff_percentage = np.sum(different_pixels) / different_pixels.size

                self.assertLessEqual(
                    diff_percentage,
                    threshold,
                    f"Images differ by {diff_percentage:.2%}, threshold is {threshold:.2%}",
                )

        except Exception as e:
            self.fail(f"Error comparing images: {e}")

    def test_histogram(self):
        result = histogram(
            project=self.project_name,
            request={
                "datasets": {
                    self.dataset_name: {
                        "compound": "DA",
                        "region": ["OF"],
                        "remove_outliers": {"grubbs": "calculated"},
                    },
                },
                "experiment": "agonist_antagonist",
            },
        )
        self.assertIsNotNone(result)
        self.assertTrue(hasattr(result, "data"))
        self.assertTrue(hasattr(result, "statistics"))

        expected_file = "./tests/results/histogram/DA in OF.png"
        actual_file = "./PROJECTS/tcb2/histogram/DA in OF.png"
        self.assert_image_similar(expected_file, actual_file, 0)

    def test_summary_histogram(self):
        result = summary_histogram(
            project=self.project_name,
            request={
                "datasets": {
                    self.dataset_name: {
                        "compound": ["DA", "NA"],
                        "region": ["OF", "PL"],
                        "remove_outliers": {"grubbs": "calculated"},
                    },
                },
                "experiment": "agonist_antagonist",
            },
        )
        self.assertIsNotNone(result)
        self.assertTrue(hasattr(result, "data"))
        self.assertTrue(hasattr(result, "statistics"))

        expected_file = "./tests/results/summary_histogram/DA, NA in OF, PL.png"
        actual_file = "./PROJECTS/tcb2/summary_histogram/DA, NA in OF, PL.png"
        self.assert_image_similar(expected_file, actual_file, 0)

    def test_correlogram(self):
        result = correlogram(
            project=self.project_name,
            request={
                "datasets": {
                    self.dataset_name: {
                        "region": "all",
                        "remove_outliers": {"grubbs": "calculated"},
                    },
                },
                "experiment": "agonist_antagonist",
            },
            between={"compound": [["5HT", "DA"]]},
        )
        self.assertIsNotNone(result)

        expected_file = "./tests/results/correlogram/all compounds in all.png"
        actual_file = "./PROJECTS/tcb2/correlogram/all compounds in all.png"
        self.assert_image_similar(expected_file, actual_file, 0)

    def test_correlogram_multiple_datasets(self):
        result = correlogram(
            project=self.project_name,
            request={
                "datasets": {
                    self.dataset_name: {
                        "region": "all",
                        "remove_outliers": {"grubbs": "calculated"},
                    },
                    "behavior": {},
                },
                "experiment": "agonist_antagonist",
            },
            between={"dataset": [["hplc", "behavior"]]},
        )
        self.assertIsNotNone(result)

        expected_file = (
            "./tests/results/correlogram/all compounds in all and all measures.png"
        )
        actual_file = (
            "./PROJECTS/tcb2/correlogram/all compounds in all and all measures.png"
        )
        self.assert_image_similar(expected_file, actual_file, 0)

    def test_network(self):
        result = network(
            project=self.project_name,
            request={
                "datasets": {
                    self.dataset_name: {
                        "region": "all",
                        "compound": ["5HT", "DA"],
                    },
                },
            },
            between={"compound": [["5HT", "DA"]]},
            layout="all_regions",
        )
        self.assertIsNotNone(result)

        expected_file = "./tests/results/network/5HT, DA in all.png"
        actual_file = "./PROJECTS/tcb2/network/5HT, DA in all.png"
        self.assert_image_similar(expected_file, actual_file, 0.015)

    def test_network_circular(self):
        result = network(
            project=self.project_name,
            request={
                "datasets": {
                    self.dataset_name: {
                        "region": "all",
                        "compound": ["5HT", "DA"],
                    },
                },
            },
            between={"compound": [["5HT", "DA"]]},
            filename="network_circular",
        )
        self.assertIsNotNone(result)

        expected_file = "./tests/results/network/network_circular.png"
        actual_file = "./PROJECTS/tcb2/network/network_circular.png"
        self.assert_image_similar(expected_file, actual_file, 0.05)

    def test_network_degrees(self):
        result = network_degrees(
            project=self.project_name,
            request={
                "datasets": {
                    self.dataset_name: {
                        "region": "all",
                    },
                },
                "experiment": "agonist_antagonist",
            },
            between={"compound": [["DA", "5HT"]]},
            custom_params={"width": 15, "height": 120},
        )
        self.assertIsNotNone(result)

        expected_file = "./tests/results/network_degrees/all compounds in all.png"
        actual_file = "./PROJECTS/tcb2/network_degrees/all compounds in all.png"
        self.assert_image_similar(expected_file, actual_file, 0.01)

    def test_network_summary(self):
        result = network_summary(
            project=self.project_name,
            request={
                "datasets": {
                    self.dataset_name: {
                        "remove_outliers": {"grubbs": "calculated"},
                    },
                },
            },
            filename="5HT-DA max degrees all regions",
            between={"compound": [["DA", "5HT"]]},
            measurement="max_degree",
            custom_params={"size": 15},
        )
        self.assertIsNotNone(result)

        expected_file = (
            "./tests/results/network_summary/5HT-DA max degrees all regions.png"
        )
        actual_file = (
            "./PROJECTS/tcb2/network_summary/5HT-DA max degrees all regions.png"
        )
        self.assert_image_similar(expected_file, actual_file, 0)

    def test_correlation(self):
        correlation(
            project="TCB2",
            x={
                "dataset": "hplc",
                "compound": "DA",
                "region": "SN",
                "remove_outliers": {"grubbs": "calculated"},
            },
            y={
                "dataset": "hplc",
                "compound": "5HT",
                "region": "OF",
                "remove_outliers": {"grubbs": "calculated"},
            },
            grouper={"group_name": "vehicles"},
        )
        # Note: correlation might not return a result object, adjust as needed

        expected_file = "./tests/results/correlation/vehicles.png"
        actual_file = "./PROJECTS/tcb2/correlation/vehicles.png"
        self.assert_image_similar(expected_file, actual_file, 0.1)

    def test_statistics_table(self):
        result = statistics_table(
            project=self.project_name,
            request={
                "datasets": {
                    self.dataset_name: {},
                },
                "experiment": "agonist_antagonist",
            },
        )
        self.assertIsNotNone(result)

        expected_file = (
            "./tests/results/statistics_table/all compounds in all regions.xlsx"
        )
        actual_file = (
            "./PROJECTS/tcb2/statistics_table/all compounds in all regions.xlsx"
        )
        actual_df = pd.read_excel(actual_file)
        expected_df = pd.read_excel(expected_file)
        pd.testing.assert_frame_equal(actual_df, expected_df)


if __name__ == "__main__":
    unittest.main()
