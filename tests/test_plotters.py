import os
import unittest
from unittest.mock import patch
from matplotlib.testing.compare import compare_images
import pandas as pd
from tests import _mpl_test_setup  # noqa: F401
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
        cls.project_name = "tcb2_test_project"
        cls.dataset_name = "hplc"

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

    def assert_image_similar(self, expected_path, actual_path, tol=2.0):
        """
        Assert two PNGs are visually similar using Matplotlib's standard comparator.

        Args:
            expected_path: baseline image path (PNG)
            actual_path:   test output image path (PNG)
            tol: RMS tolerance in pixel intensity (float).
                0 means identical, larger tolerates minor AA/font rasterization diffs.
                Typical stable values: 1.0–5.0 depending on plot complexity.
        """
        self.assertTrue(
            os.path.exists(expected_path), f"Expected not found: {expected_path}"
        )
        self.assertTrue(os.path.exists(actual_path), f"Actual not found: {actual_path}")

        res = compare_images(expected_path, actual_path, tol=tol)
        if res is not None:
            self.fail(res)

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
        actual_file = f"./PROJECTS/{self.project_name}/histogram/DA in OF.png"
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
        actual_file = (
            f"./PROJECTS/{self.project_name}/summary_histogram/DA, NA in OF, PL.png"
        )
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
        actual_file = (
            f"./PROJECTS/{self.project_name}/correlogram/all compounds in all.png"
        )
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
        actual_file = f"./PROJECTS/{self.project_name}/correlogram/all compounds in all and all measures.png"
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
        actual_file = f"./PROJECTS/{self.project_name}/network/5HT, DA in all.png"
        self.assert_image_similar(expected_file, actual_file, 4)

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
        actual_file = f"./PROJECTS/{self.project_name}/network/network_circular.png"
        self.assert_image_similar(expected_file, actual_file, 5)

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
        actual_file = (
            f"./PROJECTS/{self.project_name}/network_degrees/all compounds in all.png"
        )
        self.assert_image_similar(expected_file, actual_file, 10)

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
        actual_file = f"./PROJECTS/{self.project_name}/network_summary/5HT-DA max degrees all regions.png"
        self.assert_image_similar(expected_file, actual_file, 0.08)

    def test_correlation(self):
        correlation(
            project=self.project_name,
            x={
                "dataset": self.dataset_name,
                "compound": "DA",
                "region": "SN",
                "remove_outliers": {"grubbs": "calculated"},
            },
            y={
                "dataset": self.dataset_name,
                "compound": "5HT",
                "region": "OF",
                "remove_outliers": {"grubbs": "calculated"},
            },
            grouper={"group_name": "vehicles"},
        )

        expected_file = "./tests/results/correlation/vehicles.png"
        actual_file = f"./PROJECTS/{self.project_name}/correlation/vehicles.png"
        self.assert_image_similar(expected_file, actual_file, 20)

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
        actual_file = f"./PROJECTS/{self.project_name}/statistics_table/all compounds in all regions.xlsx"
        actual_df = pd.read_excel(actual_file)
        expected_df = pd.read_excel(expected_file)
        pd.testing.assert_frame_equal(actual_df, expected_df)


if __name__ == "__main__":
    unittest.main()
