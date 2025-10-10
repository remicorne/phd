import unittest
import os
import pandas as pd
from tests.utils.base import PlotterTestCase
from module.plotters import (
    histogram,
    summary_histogram,
    correlogram,
    network,
    network_summary,
    network_degrees,
    correlation,
    statistics_table,
)


BASE_EXPECTED_IMAGES_PATH = ENVIRONMENT_EXPECTED_IMAGES_PATH = (
    "./tests/functional/plotters/expected"
)

if os.getenv("GITHUB_ACTIONS") == "true":
    ENVIRONMENT_EXPECTED_IMAGES_PATH = os.path.join(
        ENVIRONMENT_EXPECTED_IMAGES_PATH, "pipeline"
    )


class TestPlotters(PlotterTestCase):
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

        expected_file = os.path.join(
            ENVIRONMENT_EXPECTED_IMAGES_PATH, "histogram/DA in OF.png"
        )
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

        expected_file = os.path.join(
            ENVIRONMENT_EXPECTED_IMAGES_PATH, "summary_histogram/DA, NA in OF, PL.png"
        )
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
                        "region": "alzehimers_regions",
                        "remove_outliers": {"grubbs": "calculated"},
                    },
                },
                "experiment": "agonist_antagonist",
            },
            between={"compound": [["DA", "5HT"]]},
        )
        self.assertIsNone(result)

        expected_file = os.path.join(
            ENVIRONMENT_EXPECTED_IMAGES_PATH,
            "correlogram/all compounds in alzehimers_regions.png",
        )
        actual_file = f"./PROJECTS/{self.project_name}/correlogram/all compounds in alzehimers_regions.png"
        self.assert_image_similar(expected_file, actual_file, 0)

    def test_correlogram_multiple_datasets(self):
        result = correlogram(
            project=self.project_name,
            request={
                "datasets": {
                    self.dataset_name: {
                        "region": "alzehimers_regions",
                        "remove_outliers": {"grubbs": "calculated"},
                    },
                    "behavior": {},
                },
                "experiment": "agonist_antagonist",
            },
            between={"dataset": [["hplc", "behavior"]]},
        )
        self.assertIsNone(result)

        expected_file = os.path.join(
            ENVIRONMENT_EXPECTED_IMAGES_PATH,
            "correlogram/all compounds in alzehimers_regions and all measures.png",
        )
        actual_file = f"./PROJECTS/{self.project_name}/correlogram/all compounds in alzehimers_regions and all measures.png"
        self.assert_image_similar(expected_file, actual_file, 0)

    def test_network(self):
        result = network(
            project=self.project_name,
            request={
                "datasets": {
                    self.dataset_name: {
                        "region": "alzehimers_regions",
                        "compound": ["5HT", "DA"],
                    },
                },
            },
            between={"compound": [["5HT", "DA"]]},
            layout="all_regions",
        )
        self.assertIsNotNone(result)

        expected_file = os.path.join(
            ENVIRONMENT_EXPECTED_IMAGES_PATH,
            "network/5HT, DA in alzehimers_regions.png",
        )
        actual_file = (
            f"./PROJECTS/{self.project_name}/network/5HT, DA in alzehimers_regions.png"
        )
        self.assert_image_similar(expected_file, actual_file, 4)

    def test_network_circular(self):
        result = network(
            project=self.project_name,
            request={
                "datasets": {
                    self.dataset_name: {
                        "region": "alzehimers_regions",
                        "compound": ["5HT", "DA"],
                    },
                },
            },
            between={"compound": [["5HT", "DA"]]},
            filename="network_circular",
        )
        self.assertIsNotNone(result)

        expected_file = os.path.join(
            ENVIRONMENT_EXPECTED_IMAGES_PATH, "network/network_circular.png"
        )
        actual_file = f"./PROJECTS/{self.project_name}/network/network_circular.png"
        self.assert_image_similar(expected_file, actual_file, 10)

    def test_network_degrees(self):
        result = network_degrees(
            project=self.project_name,
            request={
                "datasets": {
                    self.dataset_name: {
                        "region": "alzehimers_regions",
                    },
                },
                "experiment": "agonist_antagonist",
            },
            between={"compound": [["DA", "5HT"]]},
            custom_params={"fig_width": 15, "fig_height": 120},
        )
        self.assertIsNone(result)

        expected_file = os.path.join(
            ENVIRONMENT_EXPECTED_IMAGES_PATH,
            "network_degrees/all compounds in alzehimers_regions.png",
        )
        actual_file = f"./PROJECTS/{self.project_name}/network_degrees/all compounds in alzehimers_regions.png"
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
        self.assertIsNone(result)

        expected_file = os.path.join(
            ENVIRONMENT_EXPECTED_IMAGES_PATH,
            "network_summary/5HT-DA max degrees all regions.png",
        )
        actual_file = f"./PROJECTS/{self.project_name}/network_summary/5HT-DA max degrees all regions.png"
        self.assert_image_similar(expected_file, actual_file, 0.08)

    def test_correlation(self):
        correlation(
            project=self.project_name,
            x={
                "dataset": self.dataset_name,
                "compound": "5HT",
                "region": "OF",
                "group_name": "vehicles",
                "remove_outliers": {"grubbs": "calculated"},
            },
            y={
                "dataset": self.dataset_name,
                "compound": "DA",
                "region": "OF",
                "group_name": "vehicles",
                "remove_outliers": {"grubbs": "calculated"},
            },
            filename="vehicles",
        )

        expected_file = os.path.join(
            ENVIRONMENT_EXPECTED_IMAGES_PATH, "correlation/vehicles.png"
        )
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

        expected_file = os.path.join(
            BASE_EXPECTED_IMAGES_PATH,
            "statistics_table/all compounds in all regions.xlsx",
        )
        actual_file = f"./PROJECTS/{self.project_name}/statistics_table/all compounds in all regions.xlsx"
        actual_df = pd.read_excel(actual_file)
        expected_df = pd.read_excel(expected_file)
        pd.testing.assert_frame_equal(actual_df, expected_df)


if __name__ == "__main__":
    unittest.main()
